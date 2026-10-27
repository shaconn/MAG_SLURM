#!/usr/bin/env bash
#SBATCH --job-name=mag13_metapop_run
#SBATCH --output=mag13_metapop_run_%j.out
#SBATCH --error=mag13_metapop_run_%j.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 13_metapop_run.sh
# MetaPop: per-MAG depth/breadth, macrodiversity and microdiversity across
# all samples.
#
# Outputs:
#   metapop/norm.tsv
#   metapop/output/MetaPop/...
#   metapop/output/metapop.log
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "13_metapop_run"
exit_if_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"

METAPOP_DIR="${OUT_DIR}/metapop"
METAPOP_REF_DIR="${METAPOP_DIR}/reference"
BAM_DIR="${METAPOP_DIR}/bams"
NORM_DIR="${METAPOP_DIR}/norm"
NORM_FILE="${METAPOP_DIR}/norm.tsv"
METAPOP_OUT="${METAPOP_DIR}/output"

# ----------------------------------------------------------------------------
# Validate inputs
# ----------------------------------------------------------------------------
# --reference is a DIRECTORY of per-genome FASTAs. The legacy script tested it
# with `-s`, which is always false for a directory, so MetaPop never ran
# (legacy defect #1).
require_dir "${METAPOP_REF_DIR}"
REF_GENOME_COUNT=$(find "${METAPOP_REF_DIR}" -maxdepth 1 -name '*.fa' | wc -l)
[[ "${REF_GENOME_COUNT}" -gt 0 ]] || die "no reference genomes in ${METAPOP_REF_DIR}"

require_dir "${NORM_DIR}" "${BAM_DIR}"
N_NORM="$(count_files "${NORM_DIR}" '*.norm')"
[[ "${N_NORM}" -gt 0 ]] || die "no .norm files in ${NORM_DIR}; did 12_metapop_remap complete?"
find "${NORM_DIR}" -maxdepth 1 -name '*.norm' -type f | sort | xargs cat > "${NORM_FILE}"
awk -F'\t' '$2 == 0 || $2 == "" { print $1 }' "${NORM_FILE}" | while IFS= read -r s; do
  log_warn "zero reads for sample: ${s}"
done

N_BAM="$(count_files "${BAM_DIR}" '*.bam')"
[[ "${N_BAM}" -gt 0 ]] || die "no BAMs in ${BAM_DIR}"
[[ "${N_BAM}" -eq "${N_NORM}" ]] || log_warn "BAM count (${N_BAM}) != norm entries (${N_NORM})"
while IFS= read -r bam; do
  s="$(basename "${bam}" .bam)"
  awk -F'\t' -v s="${s}" '$1 == s { found = 1 } END { exit !found }' "${NORM_FILE}" \
    || die "BAM ${bam} has no entry in ${NORM_FILE}"
done < <(find "${BAM_DIR}" -maxdepth 1 -name '*.bam' -type f)

# ----------------------------------------------------------------------------
# Run MetaPop
# ----------------------------------------------------------------------------
mkdir -p "${METAPOP_OUT}"
log_info "MetaPop: ${N_BAM} samples x ${REF_GENOME_COUNT} reference genomes"
run_in_env metapop metapop \
  --input_samples "${BAM_DIR}" \
  --reference "${METAPOP_REF_DIR}" \
  --norm "${NORM_FILE}" \
  --output "${METAPOP_OUT}" \
  --threads "${THREADS}" \
  --id_min "${METAPOP_ID_MIN}" \
  --min_cov "${METAPOP_MIN_COV}" \
  --min_dep "${METAPOP_MIN_DEP}" \
  2>&1 | tee "${METAPOP_OUT}/metapop.log"

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
summary_block "13_metapop_run: summary" \
  "Samples (BAMs)" "${N_BAM}" \
  "Reference genomes" "${REF_GENOME_COUNT}" \
  "Norm file" "${NORM_FILE}" \
  "Output" "${METAPOP_OUT}" \
  "Log" "${METAPOP_OUT}/metapop.log"

mark_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"
