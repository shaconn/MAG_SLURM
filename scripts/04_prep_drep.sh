#!/usr/bin/env bash
#SBATCH --job-name=mag04_prep_drep
#SBATCH --output=mag04_prep_drep_%j.out
#SBATCH --error=mag04_prep_drep_%j.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 04_prep_drep.sh
# Stage every MetaBAT2 bin under a sample-prefixed name and build the dRep
# inputs from the per-sample CheckM2 reports.
#
# Outputs:
#   drep/all_bins/<S>__bin.<N>.fa
#   drep/genome_info.csv      genome,completeness,contamination
#   drep/genome_paths.txt     one absolute path per line (passed to dRep -g)
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "04_prep_drep" --no-conda
exit_if_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"

DREP_DIR="${OUT_DIR}/drep"
STAGING="${DREP_DIR}/all_bins"
GENOME_INFO="${DREP_DIR}/genome_info.csv"
GENOME_PATHS="${DREP_DIR}/genome_paths.txt"

# ----------------------------------------------------------------------------
# Stage bins
# ----------------------------------------------------------------------------
# Everything is written under OUT_DIR; nothing outside PROJECT_ROOT is created
# (legacy defect #9). Staging is rebuilt from scratch on every run.
rm -rf "${STAGING}"
mkdir -p "${STAGING}"
echo "genome,completeness,contamination" > "${GENOME_INFO}"
: > "${GENOME_PATHS}"

samples_with_bins=0
samples_without_bins=0
samples_missing_checkm2=0
bins_staged=0
bins_without_qc=0

while IFS= read -r SAMPLE || [[ -n "${SAMPLE}" ]]; do
  SAMPLE="${SAMPLE%$'\r'}"
  [[ -n "${SAMPLE}" ]] || continue
  BIN_DIR="${OUT_DIR}/metabat2/${SAMPLE}/bins"
  QC_REPORT="${OUT_DIR}/checkm2/${SAMPLE}/quality_report.tsv"

  n_bins="$(count_files "${BIN_DIR}" '*.fa')"
  if [[ "${n_bins}" -eq 0 ]]; then
    samples_without_bins=$(( samples_without_bins + 1 ))
    continue
  fi
  samples_with_bins=$(( samples_with_bins + 1 ))
  if [[ ! -s "${QC_REPORT}" ]]; then
    log_warn "[${SAMPLE}] ${n_bins} bins but no CheckM2 report at ${QC_REPORT}; bins excluded"
    samples_missing_checkm2=$(( samples_missing_checkm2 + 1 ))
    bins_without_qc=$(( bins_without_qc + n_bins ))
    continue
  fi

  sample_staged=0
  # One awk pass per sample: CheckM2 columns are located by header name.
  # Bins without a CheckM2 row are not staged — dRep would otherwise try to
  # run CheckM itself, which is not installed in the dRep environment.
  while IFS=$'\t' read -r bin comp cont; do
    src="${BIN_DIR}/${bin}.fa"
    if [[ ! -s "${src}" ]]; then
      log_warn "[${SAMPLE}] CheckM2 row for ${bin} but no file ${src}"
      continue
    fi
    dest="${STAGING}/${SAMPLE}__${bin}.fa"
    cp "${src}" "${dest}"
    echo "${SAMPLE}__${bin}.fa,${comp},${cont}" >> "${GENOME_INFO}"
    echo "${dest}" >> "${GENOME_PATHS}"
    sample_staged=$(( sample_staged + 1 ))
  done < <(awk -F'\t' '
    NR == 1 {
      for (i = 1; i <= NF; i++) {
        if ($i == "Name") n = i
        if ($i == "Completeness") c = i
        if ($i == "Contamination") k = i
      }
      if (!n || !c || !k) { print "bad header" > "/dev/stderr"; exit 1 }
      next
    }
    { print $n "\t" $c "\t" $k }
  ' "${QC_REPORT}")

  bins_staged=$(( bins_staged + sample_staged ))
  missing=$(( n_bins - sample_staged ))
  if [[ "${missing}" -gt 0 ]]; then
    log_warn "[${SAMPLE}] ${missing} bin(s) have no CheckM2 row; excluded"
    bins_without_qc=$(( bins_without_qc + missing ))
  fi
done < "${SAMPLE_LIST}"

genome_info_rows=$(( $(count_lines "${GENOME_INFO}") - 1 ))

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
summary_block "04_prep_drep: staging summary" \
  "Samples with bins" "${samples_with_bins}" \
  "Samples without bins" "${samples_without_bins}" \
  "Samples missing CheckM2" "${samples_missing_checkm2}" \
  "Bins excluded (no CheckM2 row)" "${bins_without_qc}" \
  "Total bins staged" "${bins_staged}" \
  "Rows in genome_info.csv" "${genome_info_rows}" \
  "Staging dir" "${STAGING}" \
  "Genome paths file" "${GENOME_PATHS}"

[[ "${bins_staged}" -gt 0 ]] || die "zero bins staged; check that 03_binning completed for at least one sample"
[[ "${genome_info_rows}" -eq "${bins_staged}" ]] || die "genome_info.csv rows (${genome_info_rows}) != staged bins (${bins_staged})"
mark_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"
