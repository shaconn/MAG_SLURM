#!/usr/bin/env bash
#SBATCH --job-name=mag05_drep
#SBATCH --output=mag05_drep_%j.out
#SBATCH --error=mag05_drep_%j.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 05_drep.sh
# dRep dereplicate on the staged bins.
#   Primary clustering at DREP_PRIMARY_ANI (Mash), secondary at
#   DREP_SECONDARY_ANI (DREP_S_ALGORITHM), quality filter -comp/-con applied
#   before clustering using the CheckM2 values in genome_info.csv.
#
# Outputs:
#   drep/output/dereplicated_genomes/*.fa
#   drep/output/data_tables/{Cdb,Wdb,Widb,genomeInfo}.csv
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "05_drep"
exit_if_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"

DREP_DIR="${OUT_DIR}/drep"
DREP_OUT="${DREP_DIR}/output"
GENOME_INFO="${DREP_DIR}/genome_info.csv"
GENOME_PATHS="${DREP_DIR}/genome_paths.txt"
REPS_DIR="$(drep_reps_dir)"

require_file "${GENOME_INFO}" "${GENOME_PATHS}"
N_INPUT="$(count_lines "${GENOME_PATHS}")"
[[ "${N_INPUT}" -gt 0 ]] || die "no genomes listed in ${GENOME_PATHS}"
log_info "dRep on ${N_INPUT} bins"

# An unfinished previous attempt (no checkpoint) is set aside so its partial
# comparison tables are not mixed into this run.
set_aside_dir "${DREP_OUT}"

# ----------------------------------------------------------------------------
# Run dRep
# ----------------------------------------------------------------------------
# Genomes are passed as a file of paths, not a glob: at 10k+ bins the glob
# exceeds ARG_MAX (legacy defect #3). dRep >= 3 treats a single non-FASTA
# argument to -g as a list of genome paths, one per line.
DREP_OPTS=(
  -g "${GENOME_PATHS}"
  --genomeInfo "${GENOME_INFO}"
  --S_algorithm "${DREP_S_ALGORITHM}"
  -pa "${DREP_PRIMARY_ANI}"
  -sa "${DREP_SECONDARY_ANI}"
  -nc "${DREP_MIN_COVERAGE}"
  -comp "${DREP_MIN_COMPLETENESS}"
  -con "${DREP_MAX_CONTAMINATION}"
  --multiround_primary_clustering
  -p "${THREADS}"
)
if is_true "${DREP_GREEDY_SECONDARY}"; then
  DREP_OPTS+=(--greedy_secondary_clustering --primary_chunksize "${DREP_PRIMARY_CHUNKSIZE}")
fi

run_in_env drep dRep dereplicate "${DREP_OUT}" "${DREP_OPTS[@]}"

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
require_dir "${REPS_DIR}"
N_REPS="$(count_files "${REPS_DIR}" '*.fa')"
REDUCTION="$(awk -v i="${N_INPUT}" -v r="${N_REPS}" 'BEGIN { printf "%.1f", (1 - r / i) * 100 }')"

summary_block "05_drep: dereplication summary" \
  "Input bins" "${N_INPUT}" \
  "Representatives" "${N_REPS}" \
  "Reduction" "${REDUCTION}%" \
  "Representatives dir" "${REPS_DIR}"

[[ "${N_REPS}" -gt 0 ]] || die "dRep produced zero representatives; see ${DREP_OUT}/log/"
mark_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"
