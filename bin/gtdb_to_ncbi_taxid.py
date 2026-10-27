#!/usr/bin/env python3
"""Map GTDB-Tk classifications to NCBI TaxIds and write a mettannotator samplesheet.

mettannotator needs ``prefix,assembly,taxid`` per genome, where taxid is an NCBI
TaxId. GTDB-Tk reports GTDB lineages, so each MAG's lineage is resolved here:

1. Species: GTDB genomes in the metadata TSVs that share the MAG's GTDB species
   vote on ``ncbi_species_taxid`` (majority wins).
2. Genus .. phylum: genomes in the MAG's GTDB taxon vote on the NCBI name at the
   same rank (from ``ncbi_taxonomy``); that name is looked up in the NCBI
   taxdump (names.dmp + nodes.dmp) at the matching rank. Skipped if no taxdump.
3. Fallback: the superkingdom TaxId, 2 (Bacteria) or 2157 (Archaea). The
   mettannotator docs state a higher-level TaxId is acceptable when the
   species is unknown.

Prefixes are the MAG basename, sanitised to mettannotator's allowed characters
and shortened to 24 characters (head + tail, so the bin number survives). Truncation collisions are detected and resolved
by numeric suffixing; every collision is recorded in the report.

Standard library only, so it runs with any python3 on the cluster.
"""

import argparse
import csv
import glob
import gzip
import os
import re
import sys
from collections import Counter, defaultdict

PREFIX_MAX_LEN = 24
BACTERIA_TAXID = "2"
ARCHAEA_TAXID = "2157"
RANKS = ["domain", "phylum", "class", "order", "family", "genus", "species"]
RANK_PREFIX = {"d": "domain", "p": "phylum", "c": "class", "o": "order",
               "f": "family", "g": "genus", "s": "species"}
# NCBI rank names used to look up GTDB-rank names in nodes.dmp.
NCBI_RANK = {"phylum": "phylum", "class": "class", "order": "order",
             "family": "family", "genus": "genus"}

csv.field_size_limit(sys.maxsize)


def log(msg):
    print(f"[gtdb_to_ncbi_taxid] {msg}", file=sys.stderr)


def open_text(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, newline="")


def parse_lineage(lineage):
    """'d__Bacteria;p__X;...' -> {rank: name}; empty ranks (e.g. 's__') omitted."""
    out = {}
    for part in lineage.split(";"):
        part = part.strip()
        if len(part) < 3 or part[1:3] != "__":
            continue
        rank = RANK_PREFIX.get(part[0])
        name = part[3:]
        if rank and name:
            out[rank] = name
    return out


def lineage_key(parsed, rank):
    """Lineage string truncated at `rank`, used as a join key."""
    idx = RANKS.index(rank)
    parts = []
    for r in RANKS[: idx + 1]:
        if r not in parsed:
            return None
        parts.append(parsed[r])
    return ";".join(parts)


def read_gtdbtk_summaries(paths):
    """Return {genome: (lineage, domain_taxid)} from GTDB-Tk summary TSVs."""
    mags = {}
    for path, domain_taxid in paths:
        if not os.path.exists(path):
            log(f"summary not found (skipped): {path}")
            continue
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                mags[row["user_genome"]] = (row.get("classification", ""), domain_taxid)
        log(f"read {path}")
    return mags


def load_ncbi_names(taxdump_dir):
    """Return {(rank, scientific_name): taxid} from an NCBI taxdump, or {}."""
    names_path = os.path.join(taxdump_dir, "names.dmp")
    nodes_path = os.path.join(taxdump_dir, "nodes.dmp")
    if not (os.path.exists(names_path) and os.path.exists(nodes_path)):
        log(f"NCBI taxdump not found in {taxdump_dir}; genus..phylum resolution disabled")
        return {}
    wanted_ranks = set(NCBI_RANK.values())
    rank_of = {}
    with open(nodes_path) as fh:
        for line in fh:
            f = line.rstrip("\t|\n").split("\t|\t")
            if len(f) > 2 and f[2] in wanted_ranks:
                rank_of[f[0]] = f[2]
    lookup = {}
    with open(names_path) as fh:
        for line in fh:
            f = line.rstrip("\t|\n").split("\t|\t")
            if len(f) >= 4 and f[3] == "scientific name" and f[0] in rank_of:
                lookup.setdefault((rank_of[f[0]], f[1]), f[0])
    log(f"loaded {len(lookup)} NCBI names from {taxdump_dir}")
    return lookup


def vote_from_metadata(metadata_paths, needed):
    """Stream GTDB metadata and collect votes for the taxa in `needed`.

    needed: {rank: set(lineage_key)}. Only rows inside those taxa are kept, so
    memory stays proportional to the MAG set, not the GTDB release.
    Returns {(rank, key): Counter}: species votes are NCBI species TaxIds,
    higher-rank votes are NCBI names at that rank.
    """
    votes = defaultdict(Counter)
    for path in metadata_paths:
        n = 0
        with open_text(path) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                n += 1
                parsed = parse_lineage(row.get("gtdb_taxonomy", ""))
                ncbi = parse_lineage(row.get("ncbi_taxonomy", ""))
                for rank, keys in needed.items():
                    key = lineage_key(parsed, rank)
                    if key is None or key not in keys:
                        continue
                    if rank == "species":
                        taxid = row.get("ncbi_species_taxid") or row.get("ncbi_taxid") or ""
                        if taxid.isdigit():
                            votes[(rank, key)][taxid] += 1
                    elif ncbi.get(rank):
                        votes[(rank, key)][ncbi[rank]] += 1
        log(f"scanned {n} genomes in {path}")
    return votes


def resolve(parsed, domain_taxid, votes, ncbi_names):
    """Return (taxid, resolved_rank, fallback_used)."""
    key = lineage_key(parsed, "species")
    if key and votes.get(("species", key)):
        return votes[("species", key)].most_common(1)[0][0], "species", False
    for rank in ["genus", "family", "order", "class", "phylum"]:
        key = lineage_key(parsed, rank)
        if not key or not votes.get((rank, key)):
            continue
        for name, _ in votes[(rank, key)].most_common():
            taxid = ncbi_names.get((NCBI_RANK[rank], name))
            if taxid:
                return taxid, rank, False
    return domain_taxid, "domain", True


def shorten(base, limit=PREFIX_MAX_LEN):
    """Keep the head (sample) and tail (bin number) of an over-long name, since
    dRep names are '<sample>__bin.<N>' and plain truncation drops the bin."""
    if len(base) <= limit:
        return base
    head = (limit - 1) // 2
    tail = limit - 1 - head
    return f"{base[:head]}_{base[-tail:]}"


def make_prefixes(names):
    """Map genome name -> unique mettannotator prefix (<=24 chars)."""
    prefixes = {}
    used = set()
    collisions = {}
    for name in sorted(names):
        base = re.sub(r"[^A-Za-z0-9._-]", "_", name)
        candidate = shorten(base)
        n = 1
        while candidate in used:
            n += 1
            suffix = f"_{n}"
            candidate = shorten(base, PREFIX_MAX_LEN - len(suffix)) + suffix
        if candidate != shorten(base):
            collisions[name] = candidate
        used.add(candidate)
        prefixes[name] = candidate
    return prefixes, collisions


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome-dir", required=True, help="dRep dereplicated_genomes directory")
    ap.add_argument("--extension", default="fa")
    ap.add_argument("--bac-summary", required=True)
    ap.add_argument("--ar-summary", required=True)
    ap.add_argument("--metadata-dir", required=True,
                    help="directory holding bac120_metadata_r*.tsv[.gz] / ar53_metadata_r*.tsv[.gz]")
    ap.add_argument("--taxdump-dir", default="",
                    help="NCBI taxdump directory (names.dmp, nodes.dmp); optional")
    ap.add_argument("--samplesheet", required=True)
    ap.add_argument("--report", required=True)
    args = ap.parse_args()

    genomes = sorted(glob.glob(os.path.join(args.genome_dir, f"*.{args.extension}")))
    if not genomes:
        sys.exit(f"ERROR: no *.{args.extension} genomes in {args.genome_dir}")
    names = {os.path.basename(g)[: -(len(args.extension) + 1)]: os.path.abspath(g) for g in genomes}

    mags = read_gtdbtk_summaries([(args.bac_summary, BACTERIA_TAXID),
                                  (args.ar_summary, ARCHAEA_TAXID)])
    if not mags:
        sys.exit("ERROR: neither GTDB-Tk summary exists; run 06_gtdbtk first")

    metadata = sorted(glob.glob(os.path.join(args.metadata_dir, "bac120_metadata_r*.tsv*"))
                      + glob.glob(os.path.join(args.metadata_dir, "ar53_metadata_r*.tsv*")))
    if not metadata:
        log(f"WARNING: no GTDB metadata TSVs in {args.metadata_dir}; every MAG will use the fallback")

    parsed_by_name = {}
    needed = defaultdict(set)
    for name in names:
        lineage = mags.get(name, ("", BACTERIA_TAXID))[0]
        parsed = parse_lineage(lineage)
        parsed_by_name[name] = parsed
        for rank in ["species", "genus", "family", "order", "class", "phylum"]:
            key = lineage_key(parsed, rank)
            if key:
                needed[rank].add(key)

    votes = vote_from_metadata(metadata, needed) if metadata else {}
    ncbi_names = load_ncbi_names(args.taxdump_dir) if args.taxdump_dir else {}
    prefixes, collisions = make_prefixes(names)

    rows = []
    n_fallback = 0
    n_unclassified = 0
    rank_counts = Counter()
    for name in sorted(names):
        lineage, domain_taxid = mags.get(name, ("", BACTERIA_TAXID))
        if name not in mags:
            # Not in either summary (e.g. filtered by GTDB-Tk for too few
            # marker genes). The domain is unknown; Bacteria is assumed.
            n_unclassified += 1
            taxid, rank, fallback = BACTERIA_TAXID, "unclassified", True
        else:
            taxid, rank, fallback = resolve(parsed_by_name[name], domain_taxid, votes, ncbi_names)
        n_fallback += int(fallback)
        rank_counts[rank] += 1
        rows.append({
            "genome": name,
            "prefix": prefixes[name],
            "assembly": names[name],
            "gtdb_lineage": lineage or "NA",
            "resolved_rank": rank,
            "taxid": taxid,
            "fallback_used": "yes" if fallback else "no",
            "prefix_collision": "yes" if name in collisions else "no",
        })

    os.makedirs(os.path.dirname(os.path.abspath(args.samplesheet)), exist_ok=True)
    with open(args.samplesheet, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["prefix", "assembly", "taxid"])
        for r in rows:
            w.writerow([r["prefix"], r["assembly"], r["taxid"]])
    with open(args.report, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    total = len(rows)
    log(f"genomes: {total}")
    for rank in ["species", "genus", "family", "order", "class", "phylum", "domain", "unclassified"]:
        if rank_counts[rank]:
            log(f"  resolved at {rank}: {rank_counts[rank]}")
    log(f"FALLBACK (superkingdom TaxId) used for {n_fallback}/{total} genomes "
        f"({100.0 * n_fallback / total:.1f}%)")
    if n_unclassified:
        log(f"WARNING: {n_unclassified} genomes absent from GTDB-Tk summaries; assigned TaxId 2")
    if collisions:
        log(f"prefix collisions after truncation to {PREFIX_MAX_LEN} chars: {len(collisions)} "
            "(suffixed; see prefix_collision column)")
    log(f"samplesheet: {args.samplesheet}")
    log(f"report: {args.report}")


if __name__ == "__main__":
    main()
