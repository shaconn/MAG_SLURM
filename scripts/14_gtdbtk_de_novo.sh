#!/usr/bin/env bash
#SBATCH --job-name=mag14_gtdbtk_de_novo
#SBATCH --output=mag14_gtdbtk_de_novo_%j.out
#SBATCH --error=mag14_gtdbtk_de_novo_%j.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 14_gtdbtk_de_novo.sh
# GTDB-Tk de_novo_wf on the dRep representatives, run once for bacteria and
# once for archaea into the same directory. This is the input iPHoP's
# "add_to_db" step expects for adding MAGs to its host database:
#
#   iphop add_to_db --fna_dir <OUT_DIR>/drep/output/dereplicated_genomes \
#     --gtdb_dir <OUT_DIR>/gtdbtk/iphop_de_novo \
#     --out_dir <new_db> --db_dir /full/path/to/iphop_db/<release>/
#
# Taxonomy for the rest of the pipeline still comes from 06_gtdbtk.sh
# (classify_wf); this module only builds the de novo trees. It is much slower
# than classification because the trees include all GTDB reference genomes.
#
# Outputs:
#   gtdbtk/iphop_de_novo/     de_novo_wf output (bac120 and ar53 trees)
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "14_gtdbtk_de_novo"
exit_if_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"

REPS_DIR="$(drep_reps_dir)"
DE_NOVO_OUT="${OUT_DIR}/gtdbtk/iphop_de_novo"

require_dir "${REPS_DIR}"
N_REPS="$(count_files "${REPS_DIR}" '*.fa')"
[[ "${N_REPS}" -gt 0 ]] || die "no representative genomes in ${REPS_DIR}"

GTDBTK_DATA_PATH="$(resolve_gtdbtk_data)"
export GTDBTK_DATA_PATH
mkdir -p "$(dirname "${DE_NOVO_OUT}")"
# Both runs share one directory (iPHoP reads both trees from it), so a
# partial previous attempt is set aside before starting.
set_aside_dir "${DE_NOVO_OUT}"
mkdir -p "${DE_NOVO_OUT}"
setup_workdir

# ----------------------------------------------------------------------------
# Bacterial de novo tree
# ----------------------------------------------------------------------------
log_info "GTDB-Tk de_novo_wf --bacteria on ${N_REPS} genomes (outgroup ${IPHOP_BAC_OUTGROUP})"
run_in_env gtdbtk gtdbtk de_novo_wf \
  --genome_dir "${REPS_DIR}" \
  --bacteria \
  --outgroup_taxon "${IPHOP_BAC_OUTGROUP}" \
  --out_dir "${DE_NOVO_OUT}" \
  --extension fa \
  --cpus "${THREADS}" \
  --tmpdir "${WORKDIR}" \
  --force

# ----------------------------------------------------------------------------
# Archaeal de novo tree
# ----------------------------------------------------------------------------
# A MAG set with no archaea can make this run fail. That is reported as a
# warning, not a job failure, because the bacterial tree is still valid; but
# iPHoP's add_to_db expects both trees and can then fail with an IndexError
# (the iPHoP docs suggest adding archaeal genomes, e.g. its test set).
ARCHAEA_OK=true
log_info "GTDB-Tk de_novo_wf --archaea on ${N_REPS} genomes (outgroup ${IPHOP_AR_OUTGROUP})"
if ! run_in_env gtdbtk gtdbtk de_novo_wf \
  --genome_dir "${REPS_DIR}" \
  --archaea \
  --outgroup_taxon "${IPHOP_AR_OUTGROUP}" \
  --out_dir "${DE_NOVO_OUT}" \
  --extension fa \
  --cpus "${THREADS}" \
  --tmpdir "${WORKDIR}" \
  --force; then
  ARCHAEA_OK=false
  log_warn "archaeal de_novo_wf failed (often: no archaeal MAGs); iphop add_to_db may need archaeal genomes added"
fi

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
N_BAC_TREES="$(find "${DE_NOVO_OUT}" -name '*bac120*decorated*tree*' 2>/dev/null | wc -l | tr -d ' ')"
N_AR_TREES="$(find "${DE_NOVO_OUT}" -name '*ar53*decorated*tree*' 2>/dev/null | wc -l | tr -d ' ')"
if [[ "${N_BAC_TREES}" -eq 0 ]]; then
  log_warn "no bac120 decorated tree found under ${DE_NOVO_OUT}; check ${DE_NOVO_OUT}/gtdbtk.log"
fi

summary_block "14_gtdbtk_de_novo: summary" \
  "Input genomes" "${N_REPS}" \
  "Bacterial decorated trees" "${N_BAC_TREES}" \
  "Archaeal run succeeded" "${ARCHAEA_OK}" \
  "Archaeal decorated trees" "${N_AR_TREES}" \
  "iPHoP --gtdb_dir" "${DE_NOVO_OUT}" \
  "iPHoP --fna_dir" "${REPS_DIR}"

mark_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"
