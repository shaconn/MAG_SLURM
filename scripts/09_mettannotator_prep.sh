#!/usr/bin/env bash
#SBATCH --job-name=mag09_mettannotator_prep
#SBATCH --output=mag09_mettannotator_prep_%j.out
#SBATCH --error=mag09_mettannotator_prep_%j.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 09_mettannotator_prep.sh
# Build the mettannotator samplesheet (prefix,assembly,taxid). GTDB-Tk emits
# GTDB lineages, so bin/gtdb_to_ncbi_taxid.py maps each one to an NCBI TaxId.
#
# Outputs:
#   mettannotator/samplesheet.csv
#   mettannotator/taxid_resolution_report.tsv
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "09_mettannotator_prep" --no-conda
exit_if_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"

REPS_DIR="$(drep_reps_dir)"
METT_DIR="${OUT_DIR}/mettannotator"
SAMPLESHEET="${METT_DIR}/samplesheet.csv"
REPORT="${METT_DIR}/taxid_resolution_report.tsv"
BAC="$(gtdbtk_summary bac120)"
AR="$(gtdbtk_summary ar53)"

require_dir "${REPS_DIR}"
require_cmd python3
[[ -s "${BAC}" || -s "${AR}" ]] || die "no GTDB-Tk summary found (${BAC}, ${AR}); run 06_gtdbtk first"
if ! compgen -G "${DB_ROOT}/gtdbtk/bac120_metadata_r*.tsv*" > /dev/null; then
  log_warn "GTDB metadata TSVs not found in ${DB_ROOT}/gtdbtk; all MAGs will use the superkingdom fallback (setup/download_databases.sh --only gtdbtk)"
fi
mkdir -p "${METT_DIR}"

# ----------------------------------------------------------------------------
# Map GTDB lineages to NCBI TaxIds
# ----------------------------------------------------------------------------
python3 "${MAG_SLURM_ROOT}/bin/gtdb_to_ncbi_taxid.py" \
  --genome-dir "${REPS_DIR}" \
  --extension fa \
  --bac-summary "${BAC}" \
  --ar-summary "${AR}" \
  --metadata-dir "${DB_ROOT}/gtdbtk" \
  --taxdump-dir "${DB_ROOT}/gtdbtk/ncbi_taxdump" \
  --samplesheet "${SAMPLESHEET}" \
  --report "${REPORT}"

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
require_file "${SAMPLESHEET}" "${REPORT}"
N_ROWS=$(( $(count_lines "${SAMPLESHEET}") - 1 ))
N_FALLBACK="$(awk -F'\t' 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "fallback_used") c = i; next } $c == "yes" { n++ } END { print n + 0 }' "${REPORT}")"
N_DUP="$(tail -n +2 "${SAMPLESHEET}" | cut -d, -f1 | sort | uniq -d | wc -l | tr -d ' ')"
[[ "${N_DUP}" -eq 0 ]] || die "samplesheet has ${N_DUP} duplicate prefixes; mettannotator would abort"

summary_block "09_mettannotator_prep: summary" \
  "Genomes in samplesheet" "${N_ROWS}" \
  "Superkingdom TaxId fallback" "${N_FALLBACK}" \
  "Samplesheet" "${SAMPLESHEET}" \
  "Resolution report" "${REPORT}"

mark_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"
