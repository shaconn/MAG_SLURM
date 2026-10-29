#!/usr/bin/env bash
#SBATCH --job-name=mag_build_host_indexes
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# setup/build_host_indexes.sh
# Download host references and build bowtie2 indexes under DB_ROOT/hosts/.
#
#   human/        GRCh38 no-alt analysis set
#   mouse/        GRCm39
#   phix/         phiX174 (NC_001422.1)
#   human_phix/   GRCh38 + phiX, one index  (when HOST_INCLUDE_PHIX=true)
#   mouse_phix/   GRCm39 + phiX, one index  (when HOST_INCLUDE_PHIX=true)
#
# Combined references give one removal pass against one index instead of two
# sequential passes. bowtie2-build on a mammalian genome takes hours, so by
# default this script submits itself to SLURM; --local runs it in this shell.
# Idempotent: an index that already exists is skipped.
#
# Usage:
#   setup/build_host_indexes.sh [--species human,mouse,phix] [--local] [--config FILE]
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
MODULE_NAME="build_host_indexes"
install_traps

SPECIES="human,mouse,phix"
LOCAL=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --species) SPECIES="${2:?--species needs a value}"; shift 2 ;;
    --local) LOCAL=true; shift ;;
    --config) MAG_SLURM_CONFIG="${2:?--config needs a value}"; export MAG_SLURM_CONFIG; shift 2 ;;
    -h | --help) sed -n '5,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done

load_config
HOSTS="${DB_ROOT}/hosts"
mkdir -p "${HOSTS}"

if ! is_true "${LOCAL}" && [[ -z "${SLURM_JOB_ID:-}" ]]; then
  submit_self "${MAG_SLURM_ROOT}/setup/build_host_indexes.sh" build_host_indexes 16 64G 24:00:00 --species "${SPECIES}"
  exit 0
fi

activate_conda
THREADS="${SLURM_CPUS_PER_TASK:-4}"

GRCH38_DIR="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids"
GRCH38_FILE="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
GRCM39_DIR="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39"
GRCM39_FILE="GCF_000001635.27_GRCm39_genomic.fna.gz"
PHIX_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_001422.1&rettype=fasta&retmode=text"

wanted() {
  [[ ",${SPECIES}," == *",$1,"* ]]
}

index_exists() {
  [[ -f "$1.1.bt2" || -f "$1.1.bt2l" ]]
}

# fetch_ncbi_genome <dir_url> <file> <dest.fna.gz>: verified against the
# assembly directory's md5checksums.txt.
fetch_ncbi_genome() {
  local dir_url="$1"
  local file="$2"
  local dest="$3"
  if [[ -s "${dest}" ]]; then
    log_info "reference present: ${dest}"
    return 0
  fi
  fetch "${dir_url}/${file}" "${dest}.part"
  fetch "${dir_url}/md5checksums.txt" "${dest}.md5checksums.txt"
  local expected
  expected="$(awk -v f="${file}" '{ p = $2; sub(/^\.\//, "", p) } p == f { print $1; exit }' "${dest}.md5checksums.txt")"
  [[ -n "${expected}" ]] || die "no checksum for ${file} in md5checksums.txt"
  verify_md5 "${dest}.part" "${expected}"
  mv "${dest}.part" "${dest}"
  rm -f "${dest}.md5checksums.txt"
}

# build_index <name> <fasta...>
build_index() {
  local name="$1"
  shift
  local base="${HOSTS}/${name}/${name}"
  if index_exists "${base}"; then
    log_info "index exists, skipping: ${base}"
    return 0
  fi
  mkdir -p "${HOSTS}/${name}"
  local inputs
  inputs="$(IFS=,; echo "$*")"
  log_info "bowtie2-build ${name} (${THREADS} threads); this can take hours"
  run_in_env fastp bowtie2-build --threads "${THREADS}" "${inputs}" "${base}" > "${HOSTS}/${name}/bowtie2-build.log" 2>&1
  index_exists "${base}" || die "bowtie2-build finished but no index at ${base}"
  log_info "built ${base}"
}

# ----------------------------------------------------------------------------
# References
# ----------------------------------------------------------------------------
REF="${HOSTS}/references"
mkdir -p "${REF}"

# phiX is always fetched when a combined index will be built.
NEED_PHIX=false
if wanted phix || { is_true "${HOST_INCLUDE_PHIX}" && { wanted human || wanted mouse; }; }; then
  NEED_PHIX=true
fi
if is_true "${NEED_PHIX}" && [[ ! -s "${REF}/phix.fna" ]]; then
  fetch "${PHIX_URL}" "${REF}/phix.fna.part"
  grep -q '^>NC_001422' "${REF}/phix.fna.part" || die "phiX download does not look like NC_001422.1"
  mv "${REF}/phix.fna.part" "${REF}/phix.fna"
fi
if wanted human; then
  fetch_ncbi_genome "${GRCH38_DIR}" "${GRCH38_FILE}" "${REF}/GRCh38_no_alt.fna.gz"
fi
if wanted mouse; then
  fetch_ncbi_genome "${GRCM39_DIR}" "${GRCM39_FILE}" "${REF}/GRCm39.fna.gz"
fi

# ----------------------------------------------------------------------------
# Indexes
# ----------------------------------------------------------------------------
if wanted phix; then
  build_index phix "${REF}/phix.fna"
fi
if wanted human; then
  build_index human "${REF}/GRCh38_no_alt.fna.gz"
  if is_true "${HOST_INCLUDE_PHIX}"; then
    build_index human_phix "${REF}/GRCh38_no_alt.fna.gz" "${REF}/phix.fna"
  fi
fi
if wanted mouse; then
  build_index mouse "${REF}/GRCm39.fna.gz"
  if is_true "${HOST_INCLUDE_PHIX}"; then
    build_index mouse_phix "${REF}/GRCm39.fna.gz" "${REF}/phix.fna"
  fi
fi

summary_block "build_host_indexes: done" \
  "Species" "${SPECIES}" \
  "Combined phiX indexes" "${HOST_INCLUDE_PHIX}" \
  "Index root" "${HOSTS}" \
  "Reference FASTAs" "${REF} (safe to delete once indexes are built)"
