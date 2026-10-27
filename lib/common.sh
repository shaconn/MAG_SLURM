#!/usr/bin/env bash
# ----------------------------------------------------------------------------
# lib/common.sh
# Shared functions for MAG_SLURM. Sourced by every module immediately after its
# SBATCH block, and by the orchestrator, manifest builder and setup scripts.
#
# Kept compatible with bash 3.2 so the orchestrator and manifest builder also
# run on macOS login shells; the SLURM modules themselves run on the cluster.
# ----------------------------------------------------------------------------

if [[ -n "${MAG_SLURM_COMMON_LOADED:-}" ]]; then
  return 0
fi
MAG_SLURM_COMMON_LOADED=1

# Resolve the repo root from this file's own location. BASH_SOURCE is reliable
# for a sourced file; $0 is not (under SLURM it points at the spooled copy of
# the batch script, which is why the legacy `cd "$(dirname "$0")"` failed).
MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MAG_SLURM_ROOT

MODULE_NAME="${MODULE_NAME:-mag_slurm}"
WORKDIR="${WORKDIR:-}"
# Extra temporary paths (no whitespace) removed by the EXIT trap.
CLEANUP_PATHS="${CLEANUP_PATHS:-}"

# Checkpoint key used by non-array modules.
# shellcheck disable=SC2034  # used by the modules that source this file
ALL_SAMPLES_KEY="_all"

# ----------------------------------------------------------------------------
# Logging
# ----------------------------------------------------------------------------

# All log output goes to stderr so that functions which return values on
# stdout (sample_for_task, manifest_field, ...) can log safely.
_log() {
  local level="$1"
  shift
  printf '[%s] [%s] [%s] %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" "${level}" \
    "${MODULE_NAME}" "$*" >&2
}

log_info() { _log "INFO" "$@"; }
log_warn() { _log "WARN" "$@"; }
log_error() { _log "ERROR" "$@"; }

die() {
  log_error "$@"
  exit 1
}

# ERR trap handler. Reports the failing command, its location and exit code.
# With pipefail set this fires for a failure anywhere in a pipeline, which is
# exactly the case the legacy `$?`-after-pipeline checks missed.
on_err() {
  local code="$1"
  local line="$2"
  local cmd="$3"
  trap - ERR
  log_error "command failed (exit ${code}) at ${BASH_SOURCE[1]:-?}:${line}: ${cmd}"
  # Call stack, so a failure inside a helper such as run_in_env also names
  # the module line that called it.
  local i
  for (( i = 1; i < ${#FUNCNAME[@]} - 1; i++ )); do
    log_error "  called from ${BASH_SOURCE[i + 1]}:${BASH_LINENO[i]} (${FUNCNAME[i]})"
  done
  exit "${code}"
}

# EXIT trap handler. Removes WORKDIR on every exit path: success, failure,
# and SLURM timeout/cancel (SIGTERM is converted into an exit below).
on_exit() {
  local code=$?
  if [[ -n "${WORKDIR:-}" && -d "${WORKDIR}" ]]; then
    rm -rf "${WORKDIR}" || true
  fi
  local p
  for p in ${CLEANUP_PATHS}; do
    rm -rf "${p}" || true
  done
  if [[ "${code}" -eq 0 ]]; then
    log_info "exit 0"
  else
    log_error "exit ${code}"
  fi
}

on_signal() {
  local sig="$1"
  log_error "received SIG${sig} (SLURM time limit or scancel?)"
  # 128 + signal number, the conventional shell status for a signal death.
  case "${sig}" in
    INT) exit 130 ;;
    TERM) exit 143 ;;
    *) exit 1 ;;
  esac
}

install_traps() {
  set -Eeuo pipefail
  trap 'on_err "$?" "${LINENO}" "${BASH_COMMAND}"' ERR
  trap 'on_exit' EXIT
  trap 'on_signal TERM' TERM
  trap 'on_signal INT' INT
}

# ----------------------------------------------------------------------------
# Small predicates and helpers
# ----------------------------------------------------------------------------

is_true() {
  [[ "${1:-}" == "true" ]]
}

require_cmd() {
  local cmd
  for cmd in "$@"; do
    command -v "${cmd}" >/dev/null 2>&1 || die "required command not found on PATH: ${cmd}"
  done
}

require_file() {
  local f
  for f in "$@"; do
    [[ -s "${f}" ]] || die "required file missing or empty: ${f}"
  done
}

require_dir() {
  local d
  for d in "$@"; do
    [[ -d "${d}" ]] || die "required directory missing: ${d}"
  done
}

# Count non-blank lines. Robust to a missing trailing newline, which `wc -l`
# undercounts by one.
count_lines() {
  awk 'NF { n++ } END { print n + 0 }' "$1"
}

count_files() {
  local dir="$1"
  local pattern="$2"
  find "${dir}" -maxdepth 1 -type f -name "${pattern}" 2>/dev/null | wc -l | tr -d ' '
}

# Memory available to this job in whole GB. SLURM_MEM_PER_NODE is in MB.
# Outside SLURM, falls back to parsing a config value such as "256G".
mem_gb() {
  local fallback="${1:-}"
  if [[ -n "${SLURM_MEM_PER_NODE:-}" ]]; then
    echo $(( SLURM_MEM_PER_NODE / 1024 ))
    return 0
  fi
  case "${fallback}" in
    *T) echo $(( ${fallback%T} * 1024 )) ;;
    *G) echo "${fallback%G}" ;;
    *M) echo $(( ${fallback%M} / 1024 )) ;;
    *) echo 8 ;;
  esac
}

# Move an existing, unfinished output directory out of the way so tools that
# refuse to write into an existing directory (DRAM) or that would mix stale
# state into a new run (dRep) start clean. Nothing is deleted.
set_aside_dir() {
  local dir="$1"
  if [[ -e "${dir}" ]]; then
    local aside
    aside="${dir}.incomplete.$(date '+%Y%m%d_%H%M%S')"
    log_warn "moving incomplete output aside: ${dir} -> ${aside}"
    mv "${dir}" "${aside}"
  fi
}

# Print a fixed-format summary block: summary_block "Title" key val key val ...
summary_block() {
  local title="$1"
  shift
  local rule="============================================================"
  {
    echo ""
    echo "${rule}"
    echo "${title}"
    echo "${rule}"
    while [[ $# -ge 2 ]]; do
      printf '%-30s : %s\n' "$1" "$2"
      shift 2
    done
    echo "${rule}"
  } >&2
}

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------

# Sources config/pipeline.conf from the repo root. MAG_SLURM_CONFIG overrides
# the path (the orchestrator exports it to every job so --config is honoured).
load_config() {
  local conf="${MAG_SLURM_CONFIG:-${MAG_SLURM_ROOT}/config/pipeline.conf}"
  if [[ ! -f "${conf}" ]]; then
    die "config not found: ${conf} (copy config/pipeline.conf.example to config/pipeline.conf and edit it)"
  fi
  # shellcheck source-path=SCRIPTDIR source=../config/pipeline.conf.example
  source "${conf}"
  MAG_SLURM_CONFIG="${conf}"
  export MAG_SLURM_CONFIG
}

# Variables that must be set and non-empty.
_REQUIRED_VARS="PROJECT_NAME PROJECT_ROOT READ_DIR MANIFEST SAMPLE_LIST OUT_DIR
LOG_DIR CHECKPOINT_DIR QC_READ_DIR SCRATCH_ROOT DB_ROOT CONDA_SH
SLURM_ACCOUNT SLURM_PARTITION SLURM_MAIL_USER SLURM_MAIL_TYPE ARRAY_THROTTLE
CPUS_QC_ASSEMBLY MEM_QC_ASSEMBLY TIME_QC_ASSEMBLY
CPUS_BINNING MEM_BINNING TIME_BINNING
CPUS_PREP_DREP MEM_PREP_DREP TIME_PREP_DREP
CPUS_DREP MEM_DREP TIME_DREP
CPUS_GTDBTK MEM_GTDBTK TIME_GTDBTK
CPUS_DRAM MEM_DRAM TIME_DRAM
CPUS_METTANNOTATOR MEM_METTANNOTATOR TIME_METTANNOTATOR
CPUS_METAPOP_PREP MEM_METAPOP_PREP TIME_METAPOP_PREP
CPUS_METAPOP_REMAP MEM_METAPOP_REMAP TIME_METAPOP_REMAP
CPUS_METAPOP_RUN MEM_METAPOP_RUN TIME_METAPOP_RUN
RUN_FASTP RUN_HOST_REMOVAL RUN_MICROBECENSUS RUN_METAPHLAN RUN_ASSEMBLY
RUN_QUAST RUN_BINNING RUN_DREP RUN_GTDBTK RUN_METAPOP ANNOTATION_MODE
KEEP_QC_READS HOST_SPECIES HOST_INCLUDE_PHIX
FASTP_CUT_WINDOW_SIZE FASTP_CUT_MEAN_QUALITY FASTP_LENGTH_REQUIRED
FASTP_TRIM_POLY_G FASTP_DEDUP MICROBECENSUS_INPUT
METABAT_MIN_CONTIG METABAT_MAX_P METABAT_MIN_S METABAT_MIN_CV
METABAT_MIN_CLS_SIZE DREP_S_ALGORITHM DREP_PRIMARY_ANI DREP_SECONDARY_ANI
DREP_MIN_COVERAGE DREP_MIN_COMPLETENESS DREP_MAX_CONTAMINATION
DREP_GREEDY_SECONDARY DREP_PRIMARY_CHUNKSIZE
METTANNOTATOR_VERSION METTANNOTATOR_DB METTANNOTATOR_FAST
NXF_SINGULARITY_CACHEDIR METAPOP_ID_MIN METAPOP_MIN_COV METAPOP_MIN_DEP
METAPOP_DELETE_SAM"

# Variables that must be defined but may be empty ("" = skip module load).
_DEFINED_VARS="ENV_PREFIX HOST_INDEX_CUSTOM APPTAINER_MODULE NEXTFLOW_MODULE"

_BOOLEAN_VARS="RUN_FASTP RUN_HOST_REMOVAL RUN_MICROBECENSUS RUN_METAPHLAN
RUN_ASSEMBLY RUN_QUAST RUN_BINNING RUN_DREP RUN_GTDBTK RUN_METAPOP
KEEP_QC_READS HOST_INCLUDE_PHIX FASTP_TRIM_POLY_G FASTP_DEDUP
DREP_GREEDY_SECONDARY METTANNOTATOR_FAST METAPOP_DELETE_SAM"

# validate_config [--skip-manifest]
# --skip-manifest is for 01_build_manifest.sh, which creates those files.
validate_config() {
  local check_manifest=true
  if [[ "${1:-}" == "--skip-manifest" ]]; then
    check_manifest=false
  fi

  local var
  for var in ${_REQUIRED_VARS}; do
    [[ -n "${!var:-}" ]] || die "config: ${var} is unset or empty (${MAG_SLURM_CONFIG:-config})"
  done
  for var in ${_DEFINED_VARS}; do
    [[ -n "${!var+x}" ]] || die "config: ${var} must be defined (it may be \"\")"
  done
  for var in ${_BOOLEAN_VARS}; do
    case "${!var}" in
      true | false) ;;
      *) die "config: ${var} must be true or false, got '${!var}'" ;;
    esac
  done

  case "${ANNOTATION_MODE}" in
    dram | mettannotator | both | none) ;;
    *) die "config: ANNOTATION_MODE must be dram|mettannotator|both|none, got '${ANNOTATION_MODE}'" ;;
  esac
  case "${HOST_SPECIES}" in
    human | mouse | phix | none) ;;
    custom)
      [[ -n "${HOST_INDEX_CUSTOM}" ]] || die "config: HOST_SPECIES=custom requires HOST_INDEX_CUSTOM"
      ;;
    *) die "config: HOST_SPECIES must be human|mouse|phix|none|custom, got '${HOST_SPECIES}'" ;;
  esac
  case "${MICROBECENSUS_INPUT}" in
    raw | qc) ;;
    *) die "config: MICROBECENSUS_INPUT must be raw|qc, got '${MICROBECENSUS_INPUT}'" ;;
  esac
  [[ "${ARRAY_THROTTLE}" =~ ^[0-9]+$ && "${ARRAY_THROTTLE}" -ge 1 ]] \
    || die "config: ARRAY_THROTTLE must be a positive integer, got '${ARRAY_THROTTLE}'"

  [[ -d "${PROJECT_ROOT}" ]] || die "config: PROJECT_ROOT does not exist: ${PROJECT_ROOT}"
  [[ -d "${READ_DIR}" ]] || die "config: READ_DIR does not exist: ${READ_DIR}"
  [[ -d "${SCRATCH_ROOT}" ]] || die "config: SCRATCH_ROOT does not exist: ${SCRATCH_ROOT}"
  [[ -w "${SCRATCH_ROOT}" ]] || die "config: SCRATCH_ROOT is not writable: ${SCRATCH_ROOT}"
  [[ -d "${DB_ROOT}" ]] || die "config: DB_ROOT does not exist: ${DB_ROOT}"

  if is_true "${check_manifest}"; then
    [[ -s "${MANIFEST}" ]] || die "config: MANIFEST does not exist or is empty: ${MANIFEST} (run scripts/01_build_manifest.sh)"
    [[ -s "${SAMPLE_LIST}" ]] || die "config: SAMPLE_LIST does not exist or is empty: ${SAMPLE_LIST} (run scripts/01_build_manifest.sh)"
    local n_samples n_manifest
    n_samples="$(count_lines "${SAMPLE_LIST}")"
    n_manifest="$(count_lines "${MANIFEST}")"
    if [[ "${n_samples}" -ne $(( n_manifest - 1 )) ]]; then
      die "config: SAMPLE_LIST has ${n_samples} samples but MANIFEST has $(( n_manifest - 1 )) data rows; SAMPLE_LIST must be headerless and match MANIFEST (regenerate with scripts/01_build_manifest.sh)"
    fi
    local first
    first="$(head -n 1 "${SAMPLE_LIST}" | tr -d '\r')"
    [[ "${first}" != "SampleID" ]] || die "config: SAMPLE_LIST starts with a 'SampleID' header line; it must be headerless"
    local dups
    dups="$(tr -d '\r' < "${SAMPLE_LIST}" | awk 'NF' | sort | uniq -d | head -n 5 | tr '\n' ' ')"
    [[ -z "${dups}" ]] || die "config: duplicate SampleIDs in SAMPLE_LIST: ${dups}"
  fi
}

# ----------------------------------------------------------------------------
# Environment execution
# ----------------------------------------------------------------------------

_ENV_RUNNER=""
_ENV_RUNNER_FLAG=""

activate_conda() {
  [[ -f "${CONDA_SH}" ]] || die "config: CONDA_SH not found: ${CONDA_SH}"
  # conda's shell hooks reference unset variables; relax nounset while sourcing.
  set +u
  # shellcheck source=/dev/null
  source "${CONDA_SH}"
  local mamba_sh="${CONDA_SH%/conda.sh}/mamba.sh"
  if [[ -f "${mamba_sh}" ]]; then
    # shellcheck source=/dev/null
    source "${mamba_sh}"
  fi
  set -u
  _probe_env_runner
}

# Capability probe, run once. `--no-capture-output` streams tool output live
# instead of buffering it until exit, but it is not accepted by every
# mamba/conda build (legacy defect #13). Builds without the flag (mamba 2.x /
# micromamba) stream by default.
_probe_env_runner() {
  local runner
  for runner in mamba conda micromamba; do
    if command -v "${runner}" >/dev/null 2>&1; then
      _ENV_RUNNER="${runner}"
      break
    fi
  done
  [[ -n "${_ENV_RUNNER}" ]] || die "none of mamba, conda or micromamba is available after sourcing ${CONDA_SH}"
  local help_text
  help_text="$("${_ENV_RUNNER}" run --help 2>&1 || true)"
  if [[ "${help_text}" == *--no-capture-output* ]]; then
    _ENV_RUNNER_FLAG="--no-capture-output"
  else
    _ENV_RUNNER_FLAG=""
  fi
  log_info "env runner: ${_ENV_RUNNER} run ${_ENV_RUNNER_FLAG:-(no capture flag)}"
}

# run_in_env <env_name> <cmd...>
run_in_env() {
  local env_name="${ENV_PREFIX:-}$1"
  shift
  [[ -n "${_ENV_RUNNER}" ]] || _probe_env_runner
  if [[ -n "${_ENV_RUNNER_FLAG}" ]]; then
    "${_ENV_RUNNER}" run "${_ENV_RUNNER_FLAG}" -n "${env_name}" "$@"
  else
    "${_ENV_RUNNER}" run -n "${env_name}" "$@"
  fi
}

# run_pipeline_in_env <env_name> <bash script> [args...]
# Runs a multi-stage shell pipeline inside one environment. The inner shell
# sets pipefail, so a failure in any stage (e.g. bowtie2 dying mid-stream)
# fails the whole command rather than leaving a truncated-but-valid BAM.
# Arguments are passed positionally ($1, $2, ...) to avoid quoting problems.
run_pipeline_in_env() {
  local env_name="$1"
  local script="$2"
  shift 2
  run_in_env "${env_name}" bash -c "set -Eeuo pipefail; ${script}" _ "$@"
}

# ----------------------------------------------------------------------------
# Module setup
# ----------------------------------------------------------------------------

# init_module <module_name> [--no-conda]
init_module() {
  MODULE_NAME="$1"
  install_traps
  load_config
  validate_config
  mkdir -p "${LOG_DIR}/${MODULE_NAME}" "${CHECKPOINT_DIR}/${MODULE_NAME}"
  THREADS="${SLURM_CPUS_PER_TASK:-1}"
  if [[ "${2:-}" != "--no-conda" ]]; then
    activate_conda
  fi
  log_info "start: job=${SLURM_JOB_ID:-interactive} array_task=${SLURM_ARRAY_TASK_ID:-none} threads=${THREADS} host=$(hostname)"
}

# Duplicate this job's output into a per-sample log under LOG_DIR/<module>/.
start_sample_log() {
  local sample="$1"
  local sample_log="${LOG_DIR}/${MODULE_NAME}/${sample}.log"
  log_info "per-sample log: ${sample_log}"
  exec > >(tee -a "${sample_log}") 2>&1
  log_info "sample ${sample} (job=${SLURM_JOB_ID:-interactive} task=${SLURM_ARRAY_TASK_ID:-none})"
}

# Creates a private temporary directory under SCRATCH_ROOT. The EXIT trap
# installed by init_module removes it on every exit path (legacy defect #5).
setup_workdir() {
  WORKDIR="$(mktemp -d -p "${SCRATCH_ROOT}" "mag_${MODULE_NAME}_XXXXXX")"
  export WORKDIR
  log_info "workdir: ${WORKDIR}"
}

# ----------------------------------------------------------------------------
# Sample resolution
# ----------------------------------------------------------------------------

# The single place in the codebase that turns an array index into a SampleID
# (legacy defect #6). SAMPLE_LIST is headerless, so index N is line N.
sample_for_task() {
  local idx="${1:-}"
  [[ "${idx}" =~ ^[0-9]+$ && "${idx}" -ge 1 ]] \
    || die "sample_for_task: invalid array index '${idx}' (is SLURM_ARRAY_TASK_ID set?)"
  local sample
  sample="$(sed -n "${idx}p" "${SAMPLE_LIST}" | tr -d '\r' | sed -e 's/[[:space:]]*$//')"
  [[ -n "${sample}" ]] || die "sample_for_task: no sample on line ${idx} of ${SAMPLE_LIST}"
  printf '%s\n' "${sample}"
}

# manifest_field <sample_id> <column>
# Looks the sample up by SampleID, never by line number. <column> is a
# 1-based column number or a header name (e.g. Lib2_F). Prints "" for an
# empty optional column; dies if the sample is absent.
manifest_field() {
  local sample="$1"
  local column="$2"
  local value
  if ! value="$(awk -F'\t' -v s="${sample}" -v c="${column}" '
      { sub(/\r$/, "") }
      NR == 1 {
        col = 0
        if (c ~ /^[0-9]+$/) { col = c + 0 }
        else { for (i = 1; i <= NF; i++) if ($i == c) col = i }
        if (col == 0) { exit 2 }
        next
      }
      $1 == s { print $col; found = 1; exit }
      END { if (!found) exit 3 }
    ' "${MANIFEST}")"; then
    die "manifest_field: sample '${sample}' or column '${column}' not found in ${MANIFEST}"
  fi
  printf '%s\n' "${value}"
}

# library_count <sample_id> -> 1 or 2
library_count() {
  local sample="$1"
  local l1f l1r l2f l2r
  l1f="$(manifest_field "${sample}" 2)"
  l1r="$(manifest_field "${sample}" 3)"
  l2f="$(manifest_field "${sample}" 4)"
  l2r="$(manifest_field "${sample}" 5)"
  [[ -n "${l1f}" && -n "${l1r}" ]] || die "manifest: sample '${sample}' has no complete Lib1 pair"
  if [[ -z "${l2f}" && -z "${l2r}" ]]; then
    echo 1
  elif [[ -n "${l2f}" && -n "${l2r}" ]]; then
    echo 2
  else
    die "manifest: sample '${sample}' has an incomplete Lib2 pair (Lib2_F='${l2f}' Lib2_R='${l2r}')"
  fi
}

# qc_reads <sample_id> <1|2> -> comma-separated list of post-QC files for that
# mate across all libraries, as written to QC_READ_DIR by 02_qc_assembly.sh.
qc_reads() {
  local sample="$1"
  local mate="$2"
  local n i f list=""
  n="$(library_count "${sample}")"
  for (( i = 1; i <= n; i++ )); do
    f="${QC_READ_DIR}/${sample}_L${i}_R${mate}.fastq.gz"
    require_file "${f}"
    list="${list:+${list},}${f}"
  done
  printf '%s\n' "${list}"
}

# ----------------------------------------------------------------------------
# Checkpointing
# ----------------------------------------------------------------------------

checkpoint_name() {
  printf '%s/%s/%s.done\n' "${CHECKPOINT_DIR}" "$1" "$2"
}

is_done() {
  [[ -f "$(checkpoint_name "$1" "$2")" ]]
}

mark_done() {
  local cp
  cp="$(checkpoint_name "$1" "$2")"
  mkdir -p "$(dirname "${cp}")"
  printf 'completed=%s job=%s task=%s host=%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" \
    "${SLURM_JOB_ID:-interactive}" "${SLURM_ARRAY_TASK_ID:-none}" "$(hostname)" > "${cp}"
  log_info "checkpoint written: ${cp}"
}

# Exit 0 immediately if this module/sample already completed. Resubmitting an
# array is therefore safe: finished samples are skipped (legacy defect #17).
exit_if_done() {
  local module="$1"
  local sample="$2"
  if is_done "${module}" "${sample}"; then
    log_info "already complete ($(checkpoint_name "${module}" "${sample}")); skipping"
    exit 0
  fi
}

# ----------------------------------------------------------------------------
# Shared paths and database resolution
# ----------------------------------------------------------------------------

drep_reps_dir() {
  printf '%s\n' "${OUT_DIR}/drep/output/dereplicated_genomes"
}

# Bowtie2 index basename for host removal, or "" when HOST_SPECIES=none.
# phiX is folded into the human/mouse indexes (one pass, one index).
resolve_host_index() {
  local base
  case "${HOST_SPECIES}" in
    none) echo ""; return 0 ;;
    custom) base="${HOST_INDEX_CUSTOM}" ;;
    phix) base="${DB_ROOT}/hosts/phix/phix" ;;
    human | mouse)
      if is_true "${HOST_INCLUDE_PHIX}"; then
        base="${DB_ROOT}/hosts/${HOST_SPECIES}_phix/${HOST_SPECIES}_phix"
      else
        base="${DB_ROOT}/hosts/${HOST_SPECIES}/${HOST_SPECIES}"
      fi
      ;;
  esac
  if [[ ! -f "${base}.1.bt2" && ! -f "${base}.1.bt2l" ]]; then
    die "host index not found: ${base}.1.bt2[l] (run setup/build_host_indexes.sh)"
  fi
  printf '%s\n' "${base}"
}

resolve_checkm2_db() {
  local db
  db="$(find "${DB_ROOT}/checkm2" -name '*.dmnd' -type f 2>/dev/null | sort | tail -n 1)"
  [[ -n "${db}" ]] || die "CheckM2 DIAMOND database not found under ${DB_ROOT}/checkm2 (run setup/download_databases.sh --only checkm2)"
  printf '%s\n' "${db}"
}

resolve_gtdbtk_data() {
  local data="${DB_ROOT}/gtdbtk/data"
  [[ -d "${data}" ]] || die "GTDB-Tk data not found at ${data} (run setup/download_databases.sh --only gtdbtk)"
  printf '%s\n' "${data}"
}

gtdbtk_summary() {
  # $1 = bac120 | ar53
  printf '%s\n' "${OUT_DIR}/gtdbtk/classify/gtdbtk.$1.summary.tsv"
}

# ----------------------------------------------------------------------------
# Setup-script helpers
# ----------------------------------------------------------------------------

# fetch <url> <dest>: resumable download with wget or curl.
fetch() {
  local url="$1"
  local dest="$2"
  mkdir -p "$(dirname "${dest}")"
  log_info "downloading ${url}"
  if command -v wget > /dev/null 2>&1; then
    wget -q --show-progress -c -O "${dest}" "${url}" 2>&1 || wget -c -O "${dest}" "${url}"
  elif command -v curl > /dev/null 2>&1; then
    curl -fL --retry 5 -C - -o "${dest}" "${url}"
  else
    die "neither wget nor curl is available"
  fi
}

# verify_md5 <file> <expected_md5>
verify_md5() {
  local file="$1"
  local expected="$2"
  local actual
  if command -v md5sum > /dev/null 2>&1; then
    actual="$(md5sum "${file}" | awk '{ print $1 }')"
  else
    actual="$(md5 -q "${file}")"
  fi
  [[ "${actual}" == "${expected}" ]] || die "checksum mismatch for ${file}: expected ${expected}, got ${actual}"
  log_info "checksum OK: $(basename "${file}")"
}

# submit_self <script> <cpus> <mem> <time> <args...>
# Re-submits a setup script as a SLURM job (used for multi-hour builds that
# must not run on a login node). The job runs the script with --local.
submit_self() {
  local script="$1"
  local cpus="$2"
  local mem="$3"
  local time="$4"
  shift 4
  require_cmd sbatch
  local name
  name="$(basename "${script}" .sh)"
  mkdir -p "${LOG_DIR}/setup"
  local jid
  jid="$(sbatch --parsable \
    --job-name="mag_${name}" \
    --account="${SLURM_ACCOUNT}" \
    --partition="${SLURM_PARTITION}" \
    --cpus-per-task="${cpus}" --mem="${mem}" --time="${time}" \
    --nodes=1 --ntasks=1 \
    --mail-user="${SLURM_MAIL_USER}" --mail-type="${SLURM_MAIL_TYPE}" \
    --output="${LOG_DIR}/setup/${name}_%j.out" \
    --error="${LOG_DIR}/setup/${name}_%j.err" \
    --export=ALL \
    "${script}" --local "$@")"
  jid="${jid%%;*}"
  log_info "submitted ${name} as job ${jid}; logs in ${LOG_DIR}/setup/"
}

# Loads an environment module if a name is configured. "" means the tool is
# already on PATH and the load is skipped rather than treated as an error.
load_env_module() {
  local name="$1"
  [[ -n "${name}" ]] || return 0
  if ! type module >/dev/null 2>&1; then
    die "environment module '${name}' requested but the 'module' command is unavailable (set it to \"\" in the config if the tool is on PATH)"
  fi
  # Lmod/Environment Modules scripts reference unset variables.
  set +u
  if ! module load "${name}"; then
    set -u
    die "module load ${name} failed"
  fi
  set -u
  log_info "module loaded: ${name}"
}
