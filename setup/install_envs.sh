#!/usr/bin/env bash
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# setup/install_envs.sh
# Create one conda environment per setup/envs/<name>.yml, named
# "${ENV_PREFIX}<name>". Idempotent: existing environments are skipped.
#
# Usage:
#   setup/install_envs.sh [--force] [--only NAME] [--config FILE]
#
#   --force      remove and recreate environments that already exist
#   --only NAME  create only setup/envs/NAME.yml
#
# SPAdes and MicrobeCensus each get an isolated environment: both have Python
# version conflicts with the rest of the toolchain. mettannotator does not
# appear here; it runs from its own containers (see 10_mettannotator.sh).
# ----------------------------------------------------------------------------

MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
MODULE_NAME="install_envs"
install_traps

FORCE=false
ONLY=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --force) FORCE=true; shift ;;
    --only) ONLY="${2:?--only needs a value}"; shift 2 ;;
    --config) MAG_SLURM_CONFIG="${2:?--config needs a value}"; export MAG_SLURM_CONFIG; shift 2 ;;
    -h | --help) sed -n '4,18p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

load_config
activate_conda

# ----------------------------------------------------------------------------
# Create environments
# ----------------------------------------------------------------------------
ENV_LIST="$(conda env list 2>/dev/null | awk '!/^#/ && NF { print $1 }')"
env_exists() {
  printf '%s\n' "${ENV_LIST}" | grep -qx -- "$1"
}

created=0
skipped=0
for yml in "${MAG_SLURM_ROOT}"/setup/envs/*.yml; do
  base="$(basename "${yml}" .yml)"
  if [[ -n "${ONLY}" && "${base}" != "${ONLY}" ]]; then
    continue
  fi
  env_name="${ENV_PREFIX}${base}"
  if env_exists "${env_name}"; then
    if is_true "${FORCE}"; then
      log_warn "removing existing environment ${env_name} (--force)"
      "${_ENV_RUNNER}" env remove -y -n "${env_name}"
    else
      log_info "exists, skipping: ${env_name}"
      skipped=$(( skipped + 1 ))
      continue
    fi
  fi
  log_info "creating ${env_name} from ${yml}"
  "${_ENV_RUNNER}" env create -y -n "${env_name}" -f "${yml}" \
    || "${_ENV_RUNNER}" env create -n "${env_name}" -f "${yml}"
  created=$(( created + 1 ))
done

if [[ -n "${ONLY}" && $(( created + skipped )) -eq 0 ]]; then
  die "no environment file setup/envs/${ONLY}.yml"
fi

summary_block "install_envs: summary" \
  "Runner" "${_ENV_RUNNER}" \
  "Environment prefix" "${ENV_PREFIX:-(none)}" \
  "Created" "${created}" \
  "Already present (skipped)" "${skipped}"
