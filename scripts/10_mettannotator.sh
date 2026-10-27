#!/usr/bin/env bash
#SBATCH --job-name=mag10_mettannotator
#SBATCH --output=mag10_mettannotator_%j.out
#SBATCH --error=mag10_mettannotator_%j.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 10_mettannotator.sh
# Runs the mettannotator Nextflow pipeline (Gurbich et al. 2025) on the
# samplesheet from 09. This is the only module that uses containers
# (Apptainer/Singularity). Nextflow runs inside this allocation with the
# pipeline's default local executor.
#
# Outputs:
#   mettannotator/results/<prefix>/functional_annotation/merged_gff/<prefix>_annotations.gff
#   mettannotator/results/pipeline_info/
#   mettannotator/work/            Nextflow work dir (large; safe to delete after success)
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "10_mettannotator" --no-conda
exit_if_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"

METT_DIR="${OUT_DIR}/mettannotator"
SAMPLESHEET="${METT_DIR}/samplesheet.csv"
RESULTS="${METT_DIR}/results"

require_file "${SAMPLESHEET}"

# ----------------------------------------------------------------------------
# Environment
# ----------------------------------------------------------------------------
load_env_module "${NEXTFLOW_MODULE}"
load_env_module "${APPTAINER_MODULE}"
require_cmd nextflow
if ! command -v apptainer > /dev/null 2>&1 && ! command -v singularity > /dev/null 2>&1; then
  die "neither apptainer nor singularity is on PATH (check APPTAINER_MODULE)"
fi

mkdir -p "${NXF_SINGULARITY_CACHEDIR}" "${METTANNOTATOR_DB}" "${METT_DIR}"
export NXF_SINGULARITY_CACHEDIR
# Nextflow work directories are large: keep them on scratch inside the
# project, not in the submit directory.
export NXF_WORK="${METT_DIR}/work"
export NXF_HOME="${NXF_HOME:-${DB_ROOT}/nextflow_home}"

# ----------------------------------------------------------------------------
# Run mettannotator
# ----------------------------------------------------------------------------
NF_OPTS=(
  -r "${METTANNOTATOR_VERSION}"
  -profile singularity
  -resume
  --input "${SAMPLESHEET}"
  --outdir "${RESULTS}"
  --dbs "${METTANNOTATOR_DB}"
)
if is_true "${METTANNOTATOR_FAST}"; then
  NF_OPTS+=(--fast)
fi

N_GENOMES=$(( $(count_lines "${SAMPLESHEET}") - 1 ))
log_info "mettannotator ${METTANNOTATOR_VERSION} on ${N_GENOMES} genomes (fast=${METTANNOTATOR_FAST})"
# Run from METT_DIR so .nextflow/ and .nextflow.log land beside the results;
# -resume relies on that cache being in the launch directory.
cd "${METT_DIR}"
nextflow run ebi-metagenomics/mettannotator "${NF_OPTS[@]}"

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
N_GFF="$(find "${RESULTS}" -path '*/functional_annotation/merged_gff/*_annotations.gff' 2>/dev/null | wc -l | tr -d ' ')"
summary_block "10_mettannotator: summary" \
  "Genomes submitted" "${N_GENOMES}" \
  "Merged annotation GFFs" "${N_GFF}" \
  "Results" "${RESULTS}" \
  "Nextflow work dir" "${NXF_WORK}"
[[ "${N_GFF}" -eq "${N_GENOMES}" ]] || log_warn "expected ${N_GENOMES} merged GFFs, found ${N_GFF}"

mark_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"
