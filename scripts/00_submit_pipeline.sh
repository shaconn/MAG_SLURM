#!/usr/bin/env bash
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 00_submit_pipeline.sh
# Orchestrator. Run interactively from a login node; it is not a SLURM job.
# Submits every enabled module with sbatch --parsable and chains them with
# --dependency=afterok.
#
# Usage:
#   scripts/00_submit_pipeline.sh [--dry-run] [--from NN] [--force] [--config FILE]
#
#   --dry-run      print the resolved DAG and sbatch commands; submit nothing
#   --from NN      start at module NN (02-13); earlier modules are not submitted
#   --force        delete the checkpoints of every module being submitted, so
#                  completed samples/modules are recomputed
#   --config FILE  use FILE instead of config/pipeline.conf
#
# DAG:
#   02 (array) -> 03 (array) -> 04 -> 05 -+-> 06 -> 09 -> 10
#                                         +-> 07 -> 08
#                                         +-> 11 -> 12 (array) -> 13
# ----------------------------------------------------------------------------

MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"

MODULE_NAME="00_submit_pipeline"
install_traps

usage() {
  sed -n '4,23p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

# ----------------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------------
DRY_RUN=false
FROM=2
FORCE=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run) DRY_RUN=true; shift ;;
    --from)
      [[ "${2:-}" =~ ^[0-9]+$ ]] || die "--from needs a module number (02-13)"
      FROM=$(( 10#$2 ))
      [[ "${FROM}" -ge 2 && "${FROM}" -le 13 ]] || die "--from must be between 02 and 13"
      shift 2
      ;;
    --force) FORCE=true; shift ;;
    --config) MAG_SLURM_CONFIG="$(cd "$(dirname "${2:?--config needs a value}")" && pwd)/$(basename "$2")"; export MAG_SLURM_CONFIG; shift 2 ;;
    -h | --help) usage; exit 0 ;;
    *) die "unknown argument: $1 (see --help)" ;;
  esac
done

load_config
validate_config
N_SAMPLES="$(count_lines "${SAMPLE_LIST}")"
[[ "${N_SAMPLES}" -gt 0 ]] || die "SAMPLE_LIST is empty: ${SAMPLE_LIST}"

if ! is_true "${DRY_RUN}"; then
  require_cmd sbatch
fi

# ----------------------------------------------------------------------------
# Module registry
# ----------------------------------------------------------------------------
ALL_MODULES="02 03 04 05 06 07 08 09 10 11 12 13"

mod_name() {
  case "$1" in
    02) echo "02_qc_assembly" ;;
    03) echo "03_binning" ;;
    04) echo "04_prep_drep" ;;
    05) echo "05_drep" ;;
    06) echo "06_gtdbtk" ;;
    07) echo "07_dram_annotate" ;;
    08) echo "08_dram_distill" ;;
    09) echo "09_mettannotator_prep" ;;
    10) echo "10_mettannotator" ;;
    11) echo "11_metapop_prep" ;;
    12) echo "12_metapop_remap" ;;
    13) echo "13_metapop_run" ;;
  esac
}

# Config resource suffix: CPUS_<x>, MEM_<x>, TIME_<x>. 08 shares DRAM's
# resources; 09 (a short Python job) shares PREP_DREP's.
mod_resource() {
  case "$1" in
    02) echo "QC_ASSEMBLY" ;;
    03) echo "BINNING" ;;
    04 | 09) echo "PREP_DREP" ;;
    05) echo "DREP" ;;
    06) echo "GTDBTK" ;;
    07 | 08) echo "DRAM" ;;
    10) echo "METTANNOTATOR" ;;
    11) echo "METAPOP_PREP" ;;
    12) echo "METAPOP_REMAP" ;;
    13) echo "METAPOP_RUN" ;;
  esac
}

mod_is_array() {
  case "$1" in
    02 | 03 | 12) return 0 ;;
    *) return 1 ;;
  esac
}

mod_parents() {
  case "$1" in
    02) echo "" ;;
    03) echo "02" ;;
    04) echo "03" ;;
    05) echo "04" ;;
    06 | 07 | 11) echo "05" ;;
    08) echo "07" ;;
    09) echo "06" ;;
    10) echo "09" ;;
    12) echo "11" ;;
    13) echo "12" ;;
  esac
}

annotation_has() {
  [[ "${ANNOTATION_MODE}" == "both" || "${ANNOTATION_MODE}" == "$1" ]]
}

# Prints "" if the module is enabled, otherwise the reason it is skipped.
mod_skip_reason() {
  local m="$1"
  if [[ $(( 10#${m} )) -lt "${FROM}" ]]; then
    echo "before --from $(printf '%02d' "${FROM}")"
    return 0
  fi
  case "${m}" in
    02)
      if ! is_true "${RUN_FASTP}" && ! is_true "${RUN_HOST_REMOVAL}" \
        && ! is_true "${RUN_MICROBECENSUS}" && ! is_true "${RUN_METAPHLAN}" \
        && ! is_true "${RUN_ASSEMBLY}" && ! is_true "${RUN_QUAST}"; then
        echo "all QC/assembly toggles false"
      fi
      ;;
    03) is_true "${RUN_BINNING}" || echo "RUN_BINNING=false" ;;
    04 | 05) is_true "${RUN_DREP}" || echo "RUN_DREP=false" ;;
    06) is_true "${RUN_GTDBTK}" || echo "RUN_GTDBTK=false" ;;
    07 | 08) annotation_has dram || echo "ANNOTATION_MODE=${ANNOTATION_MODE}" ;;
    09 | 10) annotation_has mettannotator || echo "ANNOTATION_MODE=${ANNOTATION_MODE}" ;;
    11 | 12 | 13) is_true "${RUN_METAPOP}" || echo "RUN_METAPOP=false" ;;
  esac
}

ENABLED=""
for m in ${ALL_MODULES}; do
  if [[ -z "$(mod_skip_reason "${m}")" ]]; then
    ENABLED="${ENABLED} ${m}"
  fi
done
[[ -n "${ENABLED}" ]] || die "no modules enabled (check toggles and --from)"

is_enabled() {
  [[ " ${ENABLED} " == *" $1 "* ]]
}

# ----------------------------------------------------------------------------
# Pre-flight checks
# ----------------------------------------------------------------------------
if ! is_true "${KEEP_QC_READS}" && { is_enabled 03 || is_enabled 12; }; then
  die "KEEP_QC_READS=false but binning (03) and/or MetaPop remap (12) need the post-QC reads in ${QC_READ_DIR}"
fi
if is_enabled 09 && ! is_enabled 06; then
  if [[ ! -s "$(gtdbtk_summary bac120)" && ! -s "$(gtdbtk_summary ar53)" ]]; then
    die "mettannotator needs GTDB-Tk taxonomy, but 06_gtdbtk is not being submitted and no previous GTDB-Tk summary exists"
  fi
fi

# ----------------------------------------------------------------------------
# Submission
# ----------------------------------------------------------------------------
# Nearest submitted ancestors of a module. A skipped parent is looked through
# to its own parents, so the dependency graph closes cleanly under any
# combination of toggles and --from.
resolve_deps() {
  local m="$1"
  local p deps=""
  for p in $(mod_parents "${m}"); do
    local jid_var="JOB_${p}"
    if [[ -n "${!jid_var:-}" ]]; then
      deps="${deps} ${!jid_var}"
    else
      deps="${deps} $(resolve_deps "${p}")"
    fi
  done
  echo "${deps}" | tr ' ' '\n' | awk 'NF && !seen[$0]++' | tr '\n' ' ' | sed 's/ *$//'
}

# Every job gets MAG_SLURM_ROOT and MAG_SLURM_CONFIG in its environment so the
# module can find lib/common.sh even though SLURM runs a spooled copy.
export MAG_SLURM_ROOT MAG_SLURM_CONFIG

DAG_LINES=""
TABLE_LINES=""
ALL_JOBS=""

submit_module() {
  local m="$1"
  local name res cpus_var mem_var time_var deps dep_names cmd jid
  name="$(mod_name "${m}")"
  res="$(mod_resource "${m}")"
  cpus_var="CPUS_${res}"
  mem_var="MEM_${res}"
  time_var="TIME_${res}"
  deps="$(resolve_deps "${m}")"

  local log_sub="%j"
  local array_desc="single"
  local -a extra=()
  if mod_is_array "${m}"; then
    log_sub="%A_%a"
    # Throttled: the legacy pipeline submitted 523 unthrottled 256 GB tasks
    # (legacy defect #4).
    extra+=(--array="1-${N_SAMPLES}%${ARRAY_THROTTLE}")
    array_desc="array 1-${N_SAMPLES}%${ARRAY_THROTTLE}"
  fi
  if [[ -n "${deps}" ]]; then
    # kill-on-invalid-dep: if an upstream job fails, dependants are cancelled
    # instead of pending forever with DependencyNeverSatisfied.
    extra+=(--dependency="afterok:${deps// /:}" --kill-on-invalid-dep=yes)
  fi

  local -a cmd_arr=(sbatch --parsable
    --job-name="mag${m}_${PROJECT_NAME}"
    --account="${SLURM_ACCOUNT}"
    --partition="${SLURM_PARTITION}"
    --cpus-per-task="${!cpus_var}"
    --mem="${!mem_var}"
    --time="${!time_var}"
    --nodes=1 --ntasks=1
    --mail-user="${SLURM_MAIL_USER}"
    --mail-type="${SLURM_MAIL_TYPE}"
    --output="${LOG_DIR}/${name}/slurm_${log_sub}.out"
    --error="${LOG_DIR}/${name}/slurm_${log_sub}.err"
    --export=ALL
    ${extra[@]+"${extra[@]}"}
    "${MAG_SLURM_ROOT}/scripts/${name}.sh")

  if is_true "${DRY_RUN}"; then
    jid="DRY${m}"
    cmd="$(printf '%q ' "${cmd_arr[@]}")"
    printf '%s\n\n' "${cmd% }"
  else
    mkdir -p "${LOG_DIR}/${name}"
    jid="$("${cmd_arr[@]}")"
    jid="${jid%%;*}"
    [[ "${jid}" =~ ^[0-9]+$ ]] || die "sbatch returned an unexpected job id for ${name}: '${jid}'"
    log_info "submitted ${name}: job ${jid}"
  fi
  printf -v "JOB_${m}" '%s' "${jid}"

  dep_names=""
  local d p
  for d in ${deps}; do
    for p in ${ALL_MODULES}; do
      local v="JOB_${p}"
      if [[ "${!v:-}" == "${d}" ]]; then
        dep_names="${dep_names:+${dep_names}, }$(mod_name "${p}")"
      fi
    done
  done
  DAG_LINES="${DAG_LINES}$(printf '  %-24s %-22s after: %s' "${name}" "[${array_desc}]" "${dep_names:-(none)}")
"
  TABLE_LINES="${TABLE_LINES}$(printf '  %-24s %s' "${name}" "${jid}")
"
  ALL_JOBS="${ALL_JOBS} ${jid}"
}

if is_true "${FORCE}"; then
  for m in ${ENABLED}; do
    cp_dir="${CHECKPOINT_DIR}/$(mod_name "${m}")"
    if is_true "${DRY_RUN}"; then
      log_info "--force: would remove checkpoints in ${cp_dir}"
    elif [[ -d "${cp_dir}" ]]; then
      log_warn "--force: removing checkpoints in ${cp_dir}"
      rm -rf "${cp_dir}"
    fi
  done
fi

if is_true "${DRY_RUN}"; then
  echo "== DRY RUN: sbatch commands (nothing submitted) =="
  echo ""
fi
for m in ${ENABLED}; do
  submit_module "${m}"
done

# ----------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------
SKIPPED_LINES=""
for m in ${ALL_MODULES}; do
  reason="$(mod_skip_reason "${m}")"
  if [[ -n "${reason}" ]]; then
    SKIPPED_LINES="${SKIPPED_LINES}$(printf '  %-24s %s' "$(mod_name "${m}")" "${reason}")
"
  fi
done

ALL_JOBS="${ALL_JOBS# }"
echo "============================================================"
echo "MAG_SLURM: ${PROJECT_NAME} ($( is_true "${DRY_RUN}" && echo "DRY RUN" || echo "submitted"))"
echo "============================================================"
echo "Config          : ${MAG_SLURM_CONFIG}"
echo "Samples         : ${N_SAMPLES}"
echo "Annotation mode : ${ANNOTATION_MODE}"
echo "MetaPop         : ${RUN_METAPOP}"
echo ""
echo "Resolved DAG:"
printf '%s' "${DAG_LINES}"
if [[ -n "${SKIPPED_LINES}" ]]; then
  echo ""
  echo "Not submitted:"
  printf '%s' "${SKIPPED_LINES}"
fi
echo ""
echo "Job IDs:"
printf '%s' "${TABLE_LINES}"
echo ""
echo "Monitor : squeue -u \"\${USER}\" -j ${ALL_JOBS// /,}"
echo "Cancel  : scancel ${ALL_JOBS}"
echo "Logs    : ${LOG_DIR}/<module>/"
echo "============================================================"
