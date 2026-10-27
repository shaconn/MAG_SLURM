# docs/legacy

> **SUPERSEDED — DO NOT EXECUTE.**
>
> `original_pipeline.md` is the hand-written pipeline that MAG_SLURM replaces. It is
> kept here, unmodified, for provenance and so the new modules can be diffed against
> it. It contains known defects (each is addressed in MAG_SLURM; see the bug-fix list
> in the top-level README) and hardcoded paths that belong to one user on one
> cluster. Do not copy code out of it and do not run it.
>
> This directory is excluded from linting and CI.

Rough correspondence between the legacy steps and the new modules:

| Legacy step | MAG_SLURM module |
|---|---|
| SRA download, renaming, R manifest | out of scope; `01_build_manifest.sh` |
| Fastp, MetaSPAdes, MetaPhlAn, QUAST | `02_qc_assembly.sh` |
| 01 MAG pipeline (bowtie2 → MetaBAT2 → CheckM2) | `03_binning.sh` |
| 02 Stage bins | `04_prep_drep.sh` |
| 03 dRep | `05_drep.sh` |
| 04 GTDB-Tk | `06_gtdbtk.sh` |
| 05/06 DRAM annotate + distill (one script) | `07_dram_annotate.sh`, `08_dram_distill.sh` |
| — | `09_mettannotator_prep.sh`, `10_mettannotator.sh` (new) |
| 07 MetaPop prep | `11_metapop_prep.sh` |
| 08 MetaPop remap | `12_metapop_remap.sh` |
| 09 MetaPop run | `13_metapop_run.sh` |
