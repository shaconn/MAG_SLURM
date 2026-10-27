#!/usr/bin/env bash
#SBATCH --job-name=mag08_dram_distill
#SBATCH --output=mag08_dram_distill_%j.out
#SBATCH --error=mag08_dram_distill_%j.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 08_dram_distill.sh
# DRAM distill: cross-MAG metabolic summary from 07's annotations.
#
# Outputs:
#   dram/distillation/metabolism_summary.xlsx
#   dram/distillation/product.html
#   dram/distillation/genome_stats.tsv
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "08_dram_distill"
exit_if_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"

ANNOT_OUT="${OUT_DIR}/dram/annotation"
DISTILL_OUT="${OUT_DIR}/dram/distillation"

require_file "${ANNOT_OUT}/annotations.tsv"

# ----------------------------------------------------------------------------
# Check DRAM databases
# ----------------------------------------------------------------------------
# Distill needs the DRAM distillate sheets, which are part of the database
# setup; fail early with an actionable message rather than mid-run.
if ! DRAM_CONFIG="$(run_in_env dram DRAM-setup.py print_config 2>&1)"; then
  die "DRAM-setup.py print_config failed; is the 'dram' environment installed? (setup/install_envs.sh)"
fi
N_DBS="$(printf '%s\n' "${DRAM_CONFIG}" | awk -F': ' '/ db: / && $2 != "None" && $2 != "" { n++ } END { print n + 0 }')"
if [[ "${N_DBS}" -eq 0 ]]; then
  die "DRAM databases are not configured (DRAM-setup.py print_config lists none). Run: setup/download_databases.sh --only dram"
fi

# ----------------------------------------------------------------------------
# Distill
# ----------------------------------------------------------------------------
set_aside_dir "${DISTILL_OUT}"

RNA_OPTS=()
if [[ -s "${ANNOT_OUT}/rrnas.tsv" ]]; then
  RNA_OPTS+=(--rrna_file "${ANNOT_OUT}/rrnas.tsv")
fi
if [[ -s "${ANNOT_OUT}/trnas.tsv" ]]; then
  RNA_OPTS+=(--trna_file "${ANNOT_OUT}/trnas.tsv")
fi

log_info "DRAM distill"
run_in_env dram DRAM.py distill \
  -i "${ANNOT_OUT}/annotations.tsv" \
  -o "${DISTILL_OUT}" \
  ${RNA_OPTS[@]+"${RNA_OPTS[@]}"}

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
require_file "${DISTILL_OUT}/genome_stats.tsv"
summary_block "08_dram_distill: summary" \
  "Metabolic summary" "${DISTILL_OUT}/metabolism_summary.xlsx" \
  "Interactive heatmap" "${DISTILL_OUT}/product.html" \
  "Genome stats" "${DISTILL_OUT}/genome_stats.tsv" \
  "rRNA/tRNA included" "$( [[ "${#RNA_OPTS[@]}" -gt 0 ]] && echo yes || echo no)"

mark_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"
