# MAG_SLURM

A config-driven SLURM pipeline that turns paired-end metagenomic FASTQs into
dereplicated, taxonomically classified, functionally annotated MAGs. Each sample is
QC'd, host-filtered, profiled and assembled with metaSPAdes. fastp runs twice:
once to make the QC reads used for profiling and coverage, and once more to
deduplicate those reads for assembly only. Bins come from
MetaBAT2 and are scored with CheckM2. Bins from all samples are dereplicated with
dRep and classified with GTDB-Tk. The representatives are annotated with DRAM
and/or mettannotator, and MetaPop can optionally add population-genetics output.
Everything is plain bash plus conda environments, with numbered modules and
one sourced config file. Jobs are chained explicitly with SLURM dependencies,
and every sample gets its own log. A sibling of `MVP_SLURM`.

```
02 qc_assembly (array) ─► 03 binning (array) ─► 04 prep_drep ─► 05 drep ─┬─► 06 gtdbtk ─► 09 mettannotator_prep ─► 10 mettannotator
                                                                          ├─► 07 dram_annotate ─► 08 dram_distill
                                                                          └─► 11 metapop_prep ─► 12 metapop_remap (array) ─► 13 metapop_run
```

## Quickstart

```bash
git clone <repo-url> MAG_SLURM && cd MAG_SLURM
cp config/pipeline.conf.example config/pipeline.conf   # edit paths, account, partition
setup/install_envs.sh                                   # conda envs from setup/envs/*.yml
setup/download_databases.sh --submit                    # all DBs under DB_ROOT (or --only NAME)
setup/build_host_indexes.sh                             # submits itself; human/mouse/phiX bowtie2 indexes
scripts/01_build_manifest.sh                            # READ_DIR -> manifest.tsv + samples.txt
scripts/00_submit_pipeline.sh --dry-run                 # inspect the DAG
scripts/00_submit_pipeline.sh                           # submit
```

The config is needed from the first setup step onward, because `DB_ROOT`,
`CONDA_SH` and the SLURM account all come from it. `PROJECT_ROOT`, `READ_DIR`,
`SCRATCH_ROOT` and `DB_ROOT` must exist before you submit anything. FASTQs are the
entry point: put them (or symlinks to them) in `READ_DIR`. The default naming is
`<Sample>_L<lane>_R{1,2}.fastq.gz`. For other layouts use
`01_build_manifest.sh --pattern`, e.g. `--pattern '*_{1,2}.fastq.gz'` for
`fasterq-dump` output.

## Configuration reference

All settings live in `config/pipeline.conf`. It is copied from
`pipeline.conf.example`, which documents every value, and git ignores it.
Every job sources it.

| Variable | Meaning |
|---|---|
| `PROJECT_NAME` | Label used in SLURM job names |
| `PROJECT_ROOT`, `OUT_DIR` | Where all output goes (use scratch) |
| `READ_DIR` | Input FASTQs; manifest paths are relative to it |
| `MANIFEST`, `SAMPLE_LIST` | `manifest.tsv` (with header) and `samples.txt` (headerless) |
| `LOG_DIR`, `CHECKPOINT_DIR`, `QC_READ_DIR` | Logs, completion sentinels, post-QC reads |
| `SCRATCH_ROOT` | Parent of per-job temp dirs (removed on every exit path) |
| `DB_ROOT` | Root of all reference databases and host indexes |
| `CONDA_SH`, `ENV_PREFIX` | conda init script; optional env name prefix (e.g. `mag_`) |
| `SLURM_ACCOUNT`, `SLURM_PARTITION`, `SLURM_MAIL_*` | Passed to every `sbatch` |
| `ARRAY_THROTTLE` | Max concurrent tasks per array job |
| `CPUS_*`, `MEM_*`, `TIME_*` | Per-module resources (08 uses `*_DRAM`, 09 uses `*_PREP_DREP`) |
| `RUN_FASTP` … `RUN_QUAST` | Steps inside module 02 |
| `RUN_BINNING`, `RUN_DREP`, `RUN_GTDBTK`, `RUN_METAPOP` | Module toggles |
| `ANNOTATION_MODE` | `dram` \| `mettannotator` \| `both` \| `none` |
| `KEEP_QC_READS` | Keep post-QC reads (required by 03 and 12) |
| `HOST_SPECIES`, `HOST_INCLUDE_PHIX`, `HOST_INDEX_CUSTOM` | Host removal index selection |
| `FASTP_*` | Trimming parameters for fastp pass 1 (QC reads, never deduplicated) |
| `FASTP_DEDUP_ASSEMBLY` | `true` (default): fastp pass 2 deduplicates the QC reads for metaSPAdes only |
| `MICROBECENSUS_INPUT` | `raw` (default) or `qc` reads |
| `METABAT_*` | MetaBAT2 parameters |
| `DREP_*` | dRep ANI thresholds, quality filter, greedy clustering |
| `METTANNOTATOR_*`, `NXF_SINGULARITY_CACHEDIR`, `*_MODULE` | mettannotator revision, DBs, container cache, modules to load |
| `METAPOP_*` | MetaPop filters; `METAPOP_DELETE_SAM=false` keeps a SAM per sample |

## Modules

| # | Script | Type | Inputs | Outputs | Typical runtime* |
|---|---|---|---|---|---|
| 00 | `00_submit_pipeline.sh` | login node | config, `samples.txt` | SLURM jobs | seconds |
| 01 | `01_build_manifest.sh` | login node | `READ_DIR` | `manifest.tsv`, `samples.txt` | seconds |
| 02 | `02_qc_assembly.sh` | array | raw FASTQs | `fastp/`, `qcreads/` (+ `qcreads/dedup/`), `microbecensus/`, `metaphlan/`, `assemblies/` | 4–48 h / sample |
| 03 | `03_binning.sh` | array | assembly, QC reads | `assemblies/rehead/`, `metabat2/<S>/bins/`, `checkm2/<S>/` | 1–6 h / sample |
| 04 | `04_prep_drep.sh` | single | bins, CheckM2 reports | `drep/all_bins/`, `genome_info.csv`, `genome_paths.txt` | < 1 h |
| 05 | `05_drep.sh` | single | staged bins | `drep/output/dereplicated_genomes/` | 2–48 h |
| 06 | `06_gtdbtk.sh` | single | representatives | `gtdbtk/classify/*.summary.tsv` | 2–24 h |
| 07 | `07_dram_annotate.sh` | single | representatives | `dram/annotation/` | 1–5 days |
| 08 | `08_dram_distill.sh` | single | DRAM annotations | `dram/distillation/` | < 2 h |
| 09 | `09_mettannotator_prep.sh` | single | representatives, GTDB-Tk | `mettannotator/samplesheet.csv`, `taxid_resolution_report.tsv` | < 30 min |
| 10 | `10_mettannotator.sh` | single | samplesheet | `mettannotator/results/` | 1–5 days |
| 11 | `11_metapop_prep.sh` | single | representatives | `metapop/reference/`, `combined_mags.fasta`, `idx/` | < 2 h |
| 12 | `12_metapop_remap.sh` | array | QC reads, index | `metapop/bams/`, `metapop/norm/` | 1–4 h / sample |
| 13 | `13_metapop_run.sh` | single | BAMs, norm, reference | `metapop/norm.tsv`, `metapop/output/` | 4–24 h |

\*Depends on depth, sample count and the cluster; these figures are rough guides.

Once the run finishes, `bin/summarize_run.py "${OUT_DIR}"` writes `mag_summary.tsv`
with one row per representative. Each row has the dRep cluster, CheckM2
completeness, contamination, N50 and size, the GTDB-Tk classification, the
mettannotator CDS count and the number of samples where MetaPop detected the MAG.

## Output layout (`OUT_DIR`)

```
manifest.tsv, samples.txt
logs/<module>/slurm_<job>.{out,err}, logs/<module>/<sample>.log
checkpoints/<module>/<sample|_all>.done
fastp/                  <S>_L<i>.{json,html} (QC pass), <S>.dedup.{json,html} (dedup pass)
qcreads/                <S>_L<i>_R{1,2}.fastq.gz: trimmed + host-removed, not deduplicated
qcreads/dedup/          <S>_R{1,2}.fastq.gz: lane-merged and deduplicated; the assembly input
microbecensus/          <S>.microbecensus
metaphlan/              <S>.metaphlan
assemblies/             <S>_scaffolds.fasta, rehead/, assembly_stats/<S>/ (QUAST)
metabat2/<S>/           depth.txt, bins/bin.<N>.fa
checkm2/<S>/            quality_report.tsv
drep/                   all_bins/, genome_info.csv, genome_paths.txt, output/
gtdbtk/                 classify/gtdbtk.{bac120,ar53}.summary.tsv
dram/                   annotation/, distillation/
mettannotator/          samplesheet.csv, taxid_resolution_report.tsv, results/, work/
metapop/                reference/, combined_mags.fasta, idx/, bams/, norm/, norm.tsv, output/
```

## Troubleshooting

**How checkpoints work.** A module writes `checkpoints/<module>/<sample>.done` when it
finishes a sample. Single modules write `_all.done` instead. Each job checks for its
sentinel first and exits 0 if it exists. Resubmitting is therefore always safe:
finished work is skipped.

**Rerun from a module.** `scripts/00_submit_pipeline.sh --from 06` submits 06 and
everything downstream of it. Modules that already finished still skip themselves
through their checkpoints. To recompute them, for example after changing dRep
thresholds, add `--force`, which deletes the checkpoints of every module being
submitted:

```bash
scripts/00_submit_pipeline.sh --from 05 --force
```

**Clear one checkpoint by hand.** `rm "${CHECKPOINT_DIR}/03_binning/SampleX.done"`
recomputes SampleX on the next submission. `rm -r "${CHECKPOINT_DIR}/03_binning"`
recomputes every sample.

**An array task failed.** Find it with `sacct -j <jobid> --format=JobID,State,ExitCode`
and read `logs/<module>/<sample>.log`. The `[ERROR]` line names the failing command,
its line number and its exit code. Downstream jobs depend on the array with `afterok`
and are submitted with `--kill-on-invalid-dep=yes`, so a failed task cancels them.
Fix the cause, then rerun:

```bash
scripts/00_submit_pipeline.sh --from 02
```

Samples that already finished exit immediately; only the failed ones do real work.
To rerun individual tasks without the rest of the DAG, run
`sbatch --array=17,42 ... scripts/02_qc_assembly.sh` with the same resource
flags `--dry-run` prints.

**Timeouts.** On SIGTERM a job logs the signal, removes its scratch directory and
exits 143. Raise the module's `TIME_*` and resubmit.

**Partial DRAM/dRep output.** When an earlier attempt left an unfinished output
directory, it is renamed to `*.incomplete.<timestamp>` and never deleted. Remove
those directories once the rerun succeeds.

**Linting.** `bash -n` and `shellcheck -x` over `scripts/ lib/ setup/`. Skip
`docs/legacy/`.

## Legacy bug-fix checklist

Each defect from the prior pipeline (`docs/legacy/original_pipeline.md`) and its fix:

| # | Defect | Where fixed |
|---|---|---|
| 1 | MetaPop tested a directory with `-s` | `13_metapop_run.sh`: `require_dir` + genome count |
| 2 | No `pipefail`; `$?` after pipelines | `set -Eeuo pipefail` + `ERR` trap (`lib/common.sh`) |
| 3 | dRep glob exceeds `ARG_MAX` | `drep/genome_paths.txt` passed to `-g` |
| 4 | Unthrottled arrays | `--array=1-N%ARRAY_THROTTLE` on every array |
| 5 | `WORKDIR` removed only on success | `EXIT` trap in `setup_workdir`/`on_exit` |
| 6 | Header off-by-one | `sample_for_task`; manifest lookup by SampleID |
| 7 | Manifest written/read at different paths | All paths from config |
| 8 | DRAM annotate/distill in one script | `07_` and `08_` are separate files |
| 9 | Stray `mkdir` outside project | Removed; nothing is written outside `OUT_DIR` |
| 10 | Host-removal SAM written then deleted | `bowtie2 -S /dev/null` |
| 11 | `runMetaBat.sh` timestamped dir | `jgi_summarize_bam_contig_depths` + `metabat2` |
| 12 | Config via `cd "$(dirname "$0")"` | Resolved from `BASH_SOURCE`; `MAG_SLURM_ROOT` exported to jobs |
| 13 | `--no-capture-output` unsupported | Capability probe in `run_in_env` |
| 14 | `mktemp -p /scratch/${USER}` | `SCRATCH_ROOT` from config |
| 15 | `fastp --dedup` distorts coverage | Dedup only in a second fastp pass whose reads feed assembly alone |
| 16 | MicrobeCensus on trimmed reads | `MICROBECENSUS_INPUT=raw` default |
| 17 | No rerun safety | Checkpoint sentinels |

## Known limitations

- **Single-sample assembly and single-sample binning.** This is standard for large
  cohorts, but MetaBAT2 sees only one coverage column, which is a weak binning
  signal. Co-assembly and multi-sample binning are out of scope.
- **MetaBAT2 only.** There is no multi-binner consensus or bin refinement (e.g.
  DAS Tool).
- **Quality thresholds.** `-comp 50 -con 10` is MIMAG medium quality *minus* the
  tRNA/rRNA criteria. Do not describe the output as full MIMAG MQ.
- **metaSPAdes takes one paired-end library.** Two-lane samples are lane-merged
  before assembly.
- **mettannotator runs inside one allocation** (`CPUS_METTANNOTATOR`,
  `MEM_METTANNOTATOR`) with Nextflow's local executor.
- **DRAM taxonomy embedding.** 07 runs in parallel with GTDB-Tk, so GTDB taxonomy
  is embedded in DRAM output only when a previous GTDB-Tk result already exists.
- **TaxIds for mettannotator** are resolved from GTDB metadata (species) and the NCBI
  taxdump (genus to phylum). Otherwise they fall back to 2 (Bacteria) or 2157
  (Archaea). `mettannotator/taxid_resolution_report.tsv` records every decision.

## Citations

- **fastp**: Chen S, Zhou Y, Chen Y, Gu J. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34(17):i884–i890 (2018).
- **Bowtie 2**: Langmead B, Salzberg SL. Fast gapped-read alignment with Bowtie 2. *Nat Methods* 9:357–359 (2012).
- **SAMtools**: Danecek P, et al. Twelve years of SAMtools and BCFtools. *GigaScience* 10(2):giab008 (2021).
- **MicrobeCensus**: Nayfach S, Pollard KS. Average genome size estimation improves comparative metagenomics and sheds light on the functional ecology of the human microbiome. *Genome Biol* 16:51 (2015).
- **MetaPhlAn 4**: Blanco-Míguez A, et al. Extending and improving metagenomic taxonomic profiling with uncharacterized species using MetaPhlAn 4. *Nat Biotechnol* 41:1633–1644 (2023).
- **metaSPAdes**: Nurk S, Meleshko D, Korobeynikov A, Pevzner PA. metaSPAdes: a new versatile metagenomic assembler. *Genome Res* 27:824–834 (2017).
- **QUAST / MetaQUAST**: Mikheenko A, Saveliev V, Gurevich A. MetaQUAST: evaluation of metagenome assemblies. *Bioinformatics* 32(7):1088–1090 (2016).
- **MetaBAT 2**: Kang DD, et al. MetaBAT 2: an adaptive binning algorithm for robust and efficient genome reconstruction from metagenome assemblies. *PeerJ* 7:e7359 (2019).
- **CheckM2**: Chklovski A, Parks DH, Woodcroft BJ, Tyson GW. CheckM2: a rapid, scalable and accurate tool for assessing microbial genome quality using machine learning. *Nat Methods* 20:1203–1212 (2023).
- **dRep**: Olm MR, Brown CT, Brooks B, Banfield JF. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. *ISME J* 11:2864–2868 (2017).
- **skani**: Shaw J, Yu YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani. *Nat Methods* 20:1661–1665 (2023).
- **GTDB-Tk v2**: Chaumeil P-A, Mussig AJ, Hugenholtz P, Parks DH. GTDB-Tk v2: memory friendly classification with the genome taxonomy database. *Bioinformatics* 38(23):5315–5316 (2022).
- **GTDB**: Parks DH, et al. GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Res* 50(D1):D785–D794 (2022).
- **DRAM**: Shaffer M, et al. DRAM for distilling microbial metabolism to automate the curation of microbiome function. *Nucleic Acids Res* 48(16):8883–8900 (2020).
- **mettannotator**: Gurbich TA, Beracochea M, De Silva NH, Finn RD. mettannotator: a comprehensive and scalable Nextflow annotation pipeline for prokaryotic assemblies. *Bioinformatics* 41(2):btaf037 (2025).
- **Nextflow**: Di Tommaso P, et al. Nextflow enables reproducible computational workflows. *Nat Biotechnol* 35:316–319 (2017).
- **MetaPop**: Gregory AC, et al. MetaPop: a pipeline for macro- and microdiversity analyses and visualization of microbial and viral metagenome-derived populations. *Microbiome* 10:49 (2022).
- **MIMAG**: Bowers RM, et al. Minimum information about a single amplified genome (MISAG) and a metagenome-assembled genome (MIMAG) of bacteria and archaea. *Nat Biotechnol* 35:725–731 (2017).

## License

MIT. See `LICENSE`.
