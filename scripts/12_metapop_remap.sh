#!/usr/bin/env bash
#SBATCH --job-name=mag12_metapop_remap
#SBATCH --output=mag12_metapop_remap_%A_%a.out
#SBATCH --error=mag12_metapop_remap_%A_%a.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 12_metapop_remap.sh  (ARRAY: one task per line of SAMPLE_LIST)
# Map each sample's QC reads to the shared MAG reference from 11.
#
# BAMs go to one flat directory: MetaPop's --input_samples reads every *.bam
# in it, and the basename becomes the sample identifier in all output tables.
#
# Outputs:
#   metapop/bams/<S>.bam, <S>.bam.bai
#   metapop/norm/<S>.norm                  "<S>\t<total_reads>"
#   metapop/sams/<S>.sam                   only when METAPOP_DELETE_SAM=false
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "12_metapop_remap"

SAMPLE="$(sample_for_task "${SLURM_ARRAY_TASK_ID:-${1:-}}")"
exit_if_done "${MODULE_NAME}" "${SAMPLE}"
start_sample_log "${SAMPLE}"

METAPOP_DIR="${OUT_DIR}/metapop"
IDX="${METAPOP_DIR}/idx/combined_mags"
BAM_DIR="${METAPOP_DIR}/bams"
NORM_DIR="${METAPOP_DIR}/norm"
SAM_DIR="${METAPOP_DIR}/sams"
BAM="${BAM_DIR}/${SAMPLE}.bam"

if [[ ! -f "${IDX}.1.bt2" && ! -f "${IDX}.1.bt2l" ]]; then
  die "bowtie2 index not found: ${IDX}.1.bt2[l] (did 11_metapop_prep complete?)"
fi
READS_1="$(qc_reads "${SAMPLE}" 1)"
READS_2="$(qc_reads "${SAMPLE}" 2)"
mkdir -p "${BAM_DIR}" "${NORM_DIR}"
setup_workdir

# ----------------------------------------------------------------------------
# Map, convert, sort
# ----------------------------------------------------------------------------
# By default (METAPOP_DELETE_SAM=true) no SAM is ever written: bowtie2 streams
# straight into samtools, and pipefail catches a failure in any stage. With
# METAPOP_DELETE_SAM=false the SAM is written to metapop/sams/ and kept.
if is_true "${METAPOP_DELETE_SAM}"; then
  log_info "bowtie2 | samtools view | samtools sort"
  # shellcheck disable=SC2016  # single quotes intended: $1.. expand in the inner bash -c
  run_pipeline_in_env metapop '
    bowtie2 --threads "$1" -x "$2" -1 "$3" -2 "$4" \
      | samtools view -b - \
      | samtools sort -@ "$1" -T "$5/sort_tmp" -o "$6" -
  ' "${THREADS}" "${IDX}" "${READS_1}" "${READS_2}" "${WORKDIR}" "${BAM}.tmp"
else
  mkdir -p "${SAM_DIR}"
  SAM="${SAM_DIR}/${SAMPLE}.sam"
  log_info "bowtie2 -> ${SAM} (kept)"
  run_in_env metapop bowtie2 --threads "${THREADS}" -x "${IDX}" \
    -1 "${READS_1}" -2 "${READS_2}" -S "${SAM}"
  # shellcheck disable=SC2016  # single quotes intended: $1.. expand in the inner bash -c
  run_pipeline_in_env metapop '
    samtools view -b "$1" | samtools sort -@ "$2" -T "$3/sort_tmp" -o "$4" -
  ' "${SAM}" "${THREADS}" "${WORKDIR}" "${BAM}.tmp"
fi
mv "${BAM}.tmp" "${BAM}"
run_in_env metapop samtools index -@ "${THREADS}" "${BAM}"

# ----------------------------------------------------------------------------
# Normalisation entry
# ----------------------------------------------------------------------------
# flagstat's "in total" counts every record, i.e. both mates of every pair
# (mapped and unmapped). Under default bowtie2 settings there are no
# secondary alignments, so this is 2 x read pairs: a consistent library-size
# proxy across samples, which is all MetaPop's --norm requires.
TOTAL_READS="$(run_in_env metapop samtools flagstat -@ "${THREADS}" "${BAM}" | awk '/in total/ { print $1; exit }')"
[[ "${TOTAL_READS}" =~ ^[0-9]+$ ]] || die "could not parse total reads from samtools flagstat for ${BAM}"
if [[ "${TOTAL_READS}" -eq 0 ]]; then
  log_warn "${SAMPLE}: zero reads in BAM"
fi
printf '%s\t%s\n' "${SAMPLE}" "${TOTAL_READS}" > "${NORM_DIR}/${SAMPLE}.norm"

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
mark_done "${MODULE_NAME}" "${SAMPLE}"
summary_block "12_metapop_remap: ${SAMPLE}" \
  "BAM" "${BAM}" \
  "Total reads (flagstat)" "${TOTAL_READS}" \
  "Norm file" "${NORM_DIR}/${SAMPLE}.norm"
