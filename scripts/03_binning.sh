#!/usr/bin/env bash
#SBATCH --job-name=mag03_binning
#SBATCH --output=mag03_binning_%A_%a.out
#SBATCH --error=mag03_binning_%A_%a.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 03_binning.sh  (ARRAY: one task per line of SAMPLE_LIST)
# rehead contigs -> bowtie2 self-mapping -> depths -> MetaBAT2 -> CheckM2
#
# Outputs:
#   assemblies/rehead/<S>_scaffolds.fasta    contigs prefixed "<S>_"
#   metabat2/<S>/depth.txt
#   metabat2/<S>/bins/bin.<N>.fa
#   checkm2/<S>/quality_report.tsv           (absent when the sample has no bins)
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "03_binning"

SAMPLE="$(sample_for_task "${SLURM_ARRAY_TASK_ID:-${1:-}}")"
exit_if_done "${MODULE_NAME}" "${SAMPLE}"
start_sample_log "${SAMPLE}"

ASSEMBLY="${OUT_DIR}/assemblies/${SAMPLE}_scaffolds.fasta"
ASSEMBLY_REHEAD="${OUT_DIR}/assemblies/rehead/${SAMPLE}_scaffolds.fasta"
METABAT_DIR="${OUT_DIR}/metabat2/${SAMPLE}"
BIN_DIR="${METABAT_DIR}/bins"
CHECKM2_DIR="${OUT_DIR}/checkm2/${SAMPLE}"

require_file "${ASSEMBLY}"
READS_1="$(qc_reads "${SAMPLE}" 1)"
READS_2="$(qc_reads "${SAMPLE}" 2)"
setup_workdir

# ----------------------------------------------------------------------------
# 1. Rehead contigs
# ----------------------------------------------------------------------------
# Prefixing the SampleID makes contig names globally unique across samples.
# Idempotent: an existing non-empty reheaded file is reused.
mkdir -p "$(dirname "${ASSEMBLY_REHEAD}")"
if [[ -s "${ASSEMBLY_REHEAD}" ]]; then
  log_info "reheaded assembly exists, reusing: ${ASSEMBLY_REHEAD}"
else
  awk -v s="${SAMPLE}" '/^>/ { print ">" s "_" substr($0, 2); next } { print }' \
    "${ASSEMBLY}" > "${ASSEMBLY_REHEAD}.tmp"
  mv "${ASSEMBLY_REHEAD}.tmp" "${ASSEMBLY_REHEAD}"
  log_info "reheaded $(awk '/^>/ { n++ } END { print n + 0 }' "${ASSEMBLY_REHEAD}") contigs"
fi

# ----------------------------------------------------------------------------
# 2-3. Map the sample's own QC reads to its own assembly
# ----------------------------------------------------------------------------
IDX="${WORKDIR}/idx/${SAMPLE}"
BAM="${WORKDIR}/${SAMPLE}.sorted.bam"
mkdir -p "${WORKDIR}/idx"
log_info "bowtie2-build"
run_in_env metabat2 bowtie2-build --threads "${THREADS}" "${ASSEMBLY_REHEAD}" "${IDX}" > /dev/null

log_info "bowtie2 | samtools view | samtools sort"
# Inner pipefail: a bowtie2 failure fails the job instead of producing a
# truncated-but-valid BAM (legacy defect #2).
# shellcheck disable=SC2016  # single quotes intended: $1.. expand in the inner bash -c
run_pipeline_in_env metabat2 '
  bowtie2 --threads "$1" -x "$2" -1 "$3" -2 "$4" \
    | samtools view -b - \
    | samtools sort -@ "$1" -T "$5/sort_tmp" -o "$6" -
' "${THREADS}" "${IDX}" "${READS_1}" "${READS_2}" "${WORKDIR}" "${BAM}"
run_in_env metabat2 samtools index -@ "${THREADS}" "${BAM}"

# ----------------------------------------------------------------------------
# 4. MetaBAT2
# ----------------------------------------------------------------------------
# Explicit two-step call instead of runMetaBat.sh, whose timestamped output
# directory required a fragile `find | sort | tail -1` (legacy defect #11).
mkdir -p "${METABAT_DIR}"
rm -rf "${BIN_DIR}"
mkdir -p "${BIN_DIR}"
log_info "jgi_summarize_bam_contig_depths"
run_in_env metabat2 jgi_summarize_bam_contig_depths \
  --outputDepth "${METABAT_DIR}/depth.txt" "${BAM}"

log_info "metabat2"
run_in_env metabat2 metabat2 \
  -i "${ASSEMBLY_REHEAD}" \
  -a "${METABAT_DIR}/depth.txt" \
  -o "${BIN_DIR}/bin" \
  -t "${THREADS}" \
  --minContig "${METABAT_MIN_CONTIG}" \
  --maxP "${METABAT_MAX_P}" \
  --minS "${METABAT_MIN_S}" \
  --minCV "${METABAT_MIN_CV}" \
  --minClsSize "${METABAT_MIN_CLS_SIZE}"

BIN_COUNT="$(count_files "${BIN_DIR}" '*.fa')"
log_info "bins: ${BIN_COUNT}"

# ----------------------------------------------------------------------------
# 5. CheckM2
# ----------------------------------------------------------------------------
# Zero bins is not an error: sparse or low-depth samples legitimately yield none.
if [[ "${BIN_COUNT}" -eq 0 ]]; then
  log_warn "no bins for ${SAMPLE}; skipping CheckM2"
  rm -rf "${CHECKM2_DIR}"
else
  CHECKM2_DB="$(resolve_checkm2_db)"
  log_info "CheckM2 predict (${CHECKM2_DB})"
  mkdir -p "${CHECKM2_DIR}"
  run_in_env checkm2 checkm2 predict \
    --input "${BIN_DIR}" \
    --output-directory "${CHECKM2_DIR}" \
    --extension fa \
    --database_path "${CHECKM2_DB}" \
    --threads "${THREADS}" \
    --tmpdir "${WORKDIR}" \
    --force
  require_file "${CHECKM2_DIR}/quality_report.tsv"
fi

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
mark_done "${MODULE_NAME}" "${SAMPLE}"
summary_block "03_binning: ${SAMPLE}" \
  "Reheaded assembly" "${ASSEMBLY_REHEAD}" \
  "Bins" "${BIN_COUNT}" \
  "Bin directory" "${BIN_DIR}" \
  "CheckM2 report" "$( [[ "${BIN_COUNT}" -gt 0 ]] && echo "${CHECKM2_DIR}/quality_report.tsv" || echo "none (no bins)")"
