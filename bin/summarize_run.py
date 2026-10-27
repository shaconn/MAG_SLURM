#!/usr/bin/env python3
"""Collate per-MAG results from a finished MAG_SLURM run into one TSV.

One row per dRep representative:
  genome, drep_secondary_cluster, completeness, contamination, contig_n50,
  genome_size, gtdb_classification, mettannotator_cds, metapop_samples_detected

Every source is optional; a missing source yields "NA" in its columns.

Sources (relative to OUT_DIR):
  drep/output/dereplicated_genomes/*.fa         the representative set
  drep/output/data_tables/Cdb.csv               secondary cluster per genome
  checkm2/<sample>/quality_report.tsv           completeness, contamination, N50, size
  gtdbtk/classify/gtdbtk.{bac120,ar53}.summary.tsv
  mettannotator/taxid_resolution_report.tsv     genome -> mettannotator prefix
  mettannotator/results/<prefix>/functional_annotation/merged_gff/<prefix>_annotations.gff
  metapop/output/**                              best effort: the first abundance
                                                 table found (genomes x samples)

Usage:
  bin/summarize_run.py OUT_DIR [-o OUT_DIR/mag_summary.tsv]
"""

import argparse
import csv
import glob
import os
import sys

csv.field_size_limit(sys.maxsize)


def log(msg):
    print(f"[summarize_run] {msg}", file=sys.stderr)


def read_tsv(path, delimiter="\t"):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def drep_clusters(out_dir):
    path = os.path.join(out_dir, "drep/output/data_tables/Cdb.csv")
    if not os.path.exists(path):
        log(f"not found: {path}")
        return {}
    return {os.path.splitext(r["genome"])[0]: r.get("secondary_cluster", "NA")
            for r in read_tsv(path, ",")}


def checkm2_stats(out_dir):
    stats = {}
    for path in glob.glob(os.path.join(out_dir, "checkm2/*/quality_report.tsv")):
        sample = os.path.basename(os.path.dirname(path))
        for r in read_tsv(path):
            stats[f"{sample}__{r['Name']}"] = r
    if not stats:
        log("no CheckM2 reports found")
    return stats


def gtdb_classifications(out_dir):
    cls = {}
    for domain in ("bac120", "ar53"):
        path = os.path.join(out_dir, f"gtdbtk/classify/gtdbtk.{domain}.summary.tsv")
        if os.path.exists(path):
            for r in read_tsv(path):
                cls[r["user_genome"]] = r.get("classification", "NA")
    return cls


def mettannotator_cds(out_dir):
    report = os.path.join(out_dir, "mettannotator/taxid_resolution_report.tsv")
    if not os.path.exists(report):
        return {}
    counts = {}
    for r in read_tsv(report):
        prefix = r["prefix"]
        gff = os.path.join(out_dir, "mettannotator/results", prefix,
                           "functional_annotation/merged_gff", f"{prefix}_annotations.gff")
        if not os.path.exists(gff):
            continue
        n = 0
        with open(gff) as fh:
            for line in fh:
                if line.startswith("##FASTA"):
                    break
                f = line.split("\t", 3)
                if len(f) > 2 and f[2] == "CDS":
                    n += 1
        counts[r["genome"]] = n
    return counts


def metapop_detections(out_dir, genomes):
    """Count samples with non-zero abundance per genome, from the first MetaPop
    table whose first column holds representative genome names."""
    root = os.path.join(out_dir, "metapop/output")
    candidates = sorted(glob.glob(os.path.join(root, "**", "*bundance*.tsv"), recursive=True))
    for path in candidates:
        with open(path, newline="") as fh:
            rows = list(csv.reader(fh, delimiter="\t"))
        if len(rows) < 2:
            continue
        body = [r for r in rows[1:] if r and r[0] in genomes]
        if not body:
            continue
        counts = {}
        for r in body:
            n = 0
            for v in r[1:]:
                try:
                    n += float(v) > 0
                except ValueError:
                    pass
            counts[r[0]] = n
        log(f"MetaPop detections from {path}")
        return counts
    if os.path.isdir(root):
        log("no MetaPop genome-by-sample abundance table recognised; metapop column is NA")
    return {}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()
    out_dir = os.path.abspath(args.out_dir)
    output = args.output or os.path.join(out_dir, "mag_summary.tsv")

    reps_dir = os.path.join(out_dir, "drep/output/dereplicated_genomes")
    genomes = sorted(os.path.basename(p)[:-3] for p in glob.glob(os.path.join(reps_dir, "*.fa")))
    if not genomes:
        sys.exit(f"ERROR: no representatives in {reps_dir}")

    clusters = drep_clusters(out_dir)
    checkm = checkm2_stats(out_dir)
    gtdb = gtdb_classifications(out_dir)
    cds = mettannotator_cds(out_dir)
    detections = metapop_detections(out_dir, set(genomes))

    cols = ["genome", "drep_secondary_cluster", "completeness", "contamination",
            "contig_n50", "genome_size", "gtdb_classification", "mettannotator_cds",
            "metapop_samples_detected"]
    with open(output, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for g in genomes:
            c = checkm.get(g, {})
            w.writerow([
                g,
                clusters.get(g, "NA"),
                c.get("Completeness", "NA"),
                c.get("Contamination", "NA"),
                c.get("Contig_N50", "NA"),
                c.get("Genome_Size", "NA"),
                gtdb.get(g, "NA"),
                cds.get(g, "NA"),
                detections.get(g, "NA"),
            ])
    log(f"{len(genomes)} representatives -> {output}")


if __name__ == "__main__":
    main()
