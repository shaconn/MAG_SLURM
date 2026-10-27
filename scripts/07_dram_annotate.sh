#!/usr/bin/env bash
#SBATCH --job-name=mag07_dram_annotate
#SBATCH --output=mag07_dram_annotate_%j.out
#SBATCH --error=mag07_dram_annotate_%j.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 07_dram_annotate.sh
# DRAM annotate on the dRep representatives. Distill is a separate module
# (08_dram_distill.sh); the legacy repo submitted two jobs but kept both steps
# in one script (legacy defect #8).
#
# Outputs:
#   dram/annotation/{annotations.tsv,rrnas.tsv,trnas.tsv,genes.*}
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "07_dram_annotate"
exit_if_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"

REPS_DIR="$(drep_reps_dir)"
ANNOT_OUT="${OUT_DIR}/dram/annotation"

require_dir "${REPS_DIR}"
N_REPS="$(count_files "${REPS_DIR}" '*.fa')"
[[ "${N_REPS}" -gt 0 ]] || die "no representative genomes in ${REPS_DIR}"

# ----------------------------------------------------------------------------
# Check DRAM databases
# ----------------------------------------------------------------------------
if ! DRAM_CONFIG="$(run_in_env dram DRAM-setup.py print_config 2>&1)"; then
  die "DRAM-setup.py print_config failed; is the 'dram' environment installed? (setup/install_envs.sh)"
fi
N_DBS="$(printf '%s\n' "${DRAM_CONFIG}" | awk -F': ' '/ db: / && $2 != "None" && $2 != "" { n++ } END { print n + 0 }')"
if [[ "${N_DBS}" -eq 0 ]]; then
  die "DRAM databases are not configured (DRAM-setup.py print_config lists none). Run: setup/download_databases.sh --only dram"
fi
log_info "DRAM reports ${N_DBS} configured annotation databases"

# ----------------------------------------------------------------------------
# Annotate
# ----------------------------------------------------------------------------
# DRAM refuses to write into an existing output directory.
mkdir -p "$(dirname "${ANNOT_OUT}")"
set_aside_dir "${ANNOT_OUT}"

# GTDB-Tk taxonomy is embedded when available. In the default DAG this module
# runs alongside 06_gtdbtk, so it is only present on reruns.
TAX_OPTS=()
for domain in bac120 ar53; do
  tax_file="$(gtdbtk_summary "${domain}")"
  if [[ -s "${tax_file}" ]]; then
    TAX_OPTS+=("${tax_file}")
  fi
done
if [[ "${#TAX_OPTS[@]}" -gt 0 ]]; then
  log_info "embedding GTDB-Tk taxonomy: ${TAX_OPTS[*]}"
  TAX_OPTS=(--gtdb_taxonomy "${TAX_OPTS[@]}")
else
  log_info "no GTDB-Tk summary found; annotating without taxonomy"
fi

# The glob is deliberately quoted: DRAM expands it internally.
log_info "DRAM annotate on ${N_REPS} genomes"
run_in_env dram DRAM.py annotate \
  -i "${REPS_DIR}/*.fa" \
  -o "${ANNOT_OUT}" \
  --threads "${THREADS}" \
  --min_contig_size "${METABAT_MIN_CONTIG}" \
  ${TAX_OPTS[@]+"${TAX_OPTS[@]}"}

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
require_file "${ANNOT_OUT}/annotations.tsv"
summary_block "07_dram_annotate: summary" \
  "Genomes annotated" "${N_REPS}" \
  "Annotated genes" "$(( $(count_lines "${ANNOT_OUT}/annotations.tsv") - 1 ))" \
  "Annotations" "${ANNOT_OUT}/annotations.tsv" \
  "rRNA calls" "$( [[ -s "${ANNOT_OUT}/rrnas.tsv" ]] && echo "${ANNOT_OUT}/rrnas.tsv" || echo none)" \
  "tRNA calls" "$( [[ -s "${ANNOT_OUT}/trnas.tsv" ]] && echo "${ANNOT_OUT}/trnas.tsv" || echo none)"

mark_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"
