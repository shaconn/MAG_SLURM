#!/usr/bin/env bash
#SBATCH --job-name=mag06_gtdbtk
#SBATCH --output=mag06_gtdbtk_%j.out
#SBATCH --error=mag06_gtdbtk_%j.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 06_gtdbtk.sh
# GTDB-Tk classify_wf on the dRep representatives.
#
# Outputs:
#   gtdbtk/classify/gtdbtk.bac120.summary.tsv
#   gtdbtk/classify/gtdbtk.ar53.summary.tsv    (only if archaea are present)
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "06_gtdbtk"
exit_if_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"

REPS_DIR="$(drep_reps_dir)"
GTDBTK_OUT="${OUT_DIR}/gtdbtk"
MASH_DB="${DB_ROOT}/gtdbtk/mash_db"

require_dir "${REPS_DIR}"
N_REPS="$(count_files "${REPS_DIR}" '*.fa')"
[[ "${N_REPS}" -gt 0 ]] || die "no representative genomes in ${REPS_DIR}"

GTDBTK_DATA_PATH="$(resolve_gtdbtk_data)"
export GTDBTK_DATA_PATH
mkdir -p "${GTDBTK_OUT}" "${MASH_DB}"
setup_workdir

# ----------------------------------------------------------------------------
# Run GTDB-Tk
# ----------------------------------------------------------------------------
# --mash_db under DB_ROOT caches the reference sketch for the ANI screen
# across runs. pplacer is the memory bottleneck and does not scale past a few
# threads, so it is capped at 4 regardless of CPUS_GTDBTK.
PPLACER_CPUS=$(( THREADS < 4 ? THREADS : 4 ))
log_info "GTDB-Tk classify_wf on ${N_REPS} genomes (data: ${GTDBTK_DATA_PATH})"
run_in_env gtdbtk gtdbtk classify_wf \
  --genome_dir "${REPS_DIR}" \
  --out_dir "${GTDBTK_OUT}" \
  --extension fa \
  --cpus "${THREADS}" \
  --pplacer_cpus "${PPLACER_CPUS}" \
  --mash_db "${MASH_DB}" \
  --tmpdir "${WORKDIR}"

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
BAC="$(gtdbtk_summary bac120)"
AR="$(gtdbtk_summary ar53)"
N_BAC=0
N_AR=0
[[ -s "${BAC}" ]] && N_BAC=$(( $(count_lines "${BAC}") - 1 ))
[[ -s "${AR}" ]] && N_AR=$(( $(count_lines "${AR}") - 1 ))
[[ $(( N_BAC + N_AR )) -gt 0 ]] || die "GTDB-Tk produced no summary rows; see ${GTDBTK_OUT}/gtdbtk.log"

summary_block "06_gtdbtk: classification summary" \
  "Input genomes" "${N_REPS}" \
  "Bacterial classifications" "${N_BAC}" \
  "Archaeal classifications" "${N_AR}" \
  "Not classified" "$(( N_REPS - N_BAC - N_AR ))" \
  "bac120 summary" "${BAC}" \
  "ar53 summary" "${AR}"

mark_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"
