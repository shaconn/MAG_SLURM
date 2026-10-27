#!/usr/bin/env bash
#SBATCH --job-name=mag_download_databases
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# setup/download_databases.sh
# Download every reference database under DB_ROOT, one subdirectory each.
#
# Usage:
#   setup/download_databases.sh [--only NAME] [--skip-uniref] [--submit] [--config FILE]
#
#   --only NAME     one of: checkm2 gtdbtk dram metaphlan mettannotator
#   --skip-uniref   DRAM without UniRef90 (see note below)
#   --submit        run as a SLURM job instead of in this shell (recommended:
#                   DRAM and mettannotator take many hours and real CPU)
#
#   Database        Subdirectory        Approx. size
#   CheckM2         checkm2/            ~3 GB
#   GTDB-Tk + Mash  gtdbtk/             ~110 GB  (+ GTDB metadata, NCBI taxdump)
#   DRAM            dram/               ~500 GB
#   MetaPhlAn       metaphlan/          ~20 GB
#   mettannotator   mettannotator/      ~182 GB
#
# Environment overrides: GTDB_RELEASE (default 226), GTDB_MIRROR,
# DRAM_KEGG_LOC (path to a licensed KEGG peptide FASTA; KEGG is optional).
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
MODULE_NAME="download_databases"
install_traps

ALL_DBS="checkm2 gtdbtk dram metaphlan mettannotator"
ONLY=""
SKIP_UNIREF=false
SUBMIT=false
ORIG_ARGS=("$@")
while [[ $# -gt 0 ]]; do
  case "$1" in
    --only) ONLY="${2:?--only needs a value}"; shift 2 ;;
    --skip-uniref) SKIP_UNIREF=true; shift ;;
    --submit) SUBMIT=true; shift ;;
    --local) shift ;;
    --config) MAG_SLURM_CONFIG="${2:?--config needs a value}"; export MAG_SLURM_CONFIG; shift 2 ;;
    -h | --help) sed -n '5,27p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) die "unknown argument: $1" ;;
  esac
done
if [[ -n "${ONLY}" && " ${ALL_DBS} " != *" ${ONLY} "* ]]; then
  die "--only must be one of: ${ALL_DBS}"
fi

load_config
mkdir -p "${DB_ROOT}"

if is_true "${SUBMIT}"; then
  pass=()
  for a in "${ORIG_ARGS[@]}"; do
    if [[ "${a}" != "--submit" ]]; then
      pass+=("${a}")
    fi
  done
  submit_self "${MAG_SLURM_ROOT}/setup/download_databases.sh" 16 64G 72:00:00 ${pass[@]+"${pass[@]}"}
  exit 0
fi

activate_conda
THREADS="${SLURM_CPUS_PER_TASK:-4}"
GTDB_RELEASE="${GTDB_RELEASE:-226}"
GTDB_MIRROR="${GTDB_MIRROR:-https://data.gtdb.ecogenomic.org/releases}"

wanted() {
  [[ -z "${ONLY}" || "${ONLY}" == "$1" ]]
}

announce() {
  log_info "=== $1: $2 into $3 ==="
}

# ----------------------------------------------------------------------------
# CheckM2
# ----------------------------------------------------------------------------
if wanted checkm2; then
  dest="${DB_ROOT}/checkm2"
  announce "CheckM2 DIAMOND database" "~3 GB" "${dest}"
  if [[ -n "$(find "${dest}" -name '*.dmnd' 2>/dev/null)" ]]; then
    log_info "already present, skipping"
  else
    mkdir -p "${dest}"
    run_in_env checkm2 checkm2 database --download --path "${dest}"
  fi
  # Upstream publishes no separate checksum; a passing testrun validates the
  # database end to end.
  run_in_env checkm2 checkm2 testrun --database_path "$(resolve_checkm2_db)" --threads "${THREADS}"
fi

# ----------------------------------------------------------------------------
# GTDB-Tk reference data, Mash sketch dir, GTDB metadata, NCBI taxdump
# ----------------------------------------------------------------------------
if wanted gtdbtk; then
  dest="${DB_ROOT}/gtdbtk"
  rel_base="${GTDB_MIRROR}/release${GTDB_RELEASE}/${GTDB_RELEASE}.0"
  mkdir -p "${dest}/mash_db" "${dest}/downloads"

  # Checksums published by GTDB for this release (relative path -> md5).
  md5_file="${dest}/downloads/MD5SUM.txt"
  if ! fetch "${rel_base}/MD5SUM.txt" "${md5_file}"; then
    log_warn "could not fetch ${rel_base}/MD5SUM.txt; downloads will not be checksum-verified"
    : > "${md5_file}"
  fi
  gtdb_md5() {
    awk -v f="$1" '{ p = $2; sub(/^\.\//, "", p); sub(/^\*/, "", p) } p == f { print $1; exit }' "${md5_file}"
  }

  announce "GTDB-Tk R${GTDB_RELEASE} reference data" "~110 GB" "${dest}"
  if [[ -d "${dest}/data" ]]; then
    log_info "GTDB-Tk data already present at ${dest}/data, skipping"
  else
    pkg_rel="auxillary_files/gtdbtk_package/full_package/gtdbtk_r${GTDB_RELEASE}_data.tar.gz"
    tarball="${dest}/downloads/gtdbtk_r${GTDB_RELEASE}_data.tar.gz"
    fetch "${rel_base}/${pkg_rel}" "${tarball}"
    expected="$(gtdb_md5 "${pkg_rel}")"
    if [[ -n "${expected}" ]]; then verify_md5 "${tarball}" "${expected}"; else log_warn "no published checksum for ${pkg_rel}"; fi
    log_info "extracting (this takes a while)"
    tar -xzf "${tarball}" -C "${dest}"
    extracted="$(find "${dest}" -maxdepth 1 -mindepth 1 -type d -name 'release*' | sort | tail -n 1)"
    [[ -n "${extracted}" ]] || die "could not find the extracted release directory under ${dest}"
    ln -sfn "$(basename "${extracted}")" "${dest}/data"
    rm -f "${tarball}"
  fi

  announce "GTDB metadata (bac120, ar53)" "~1 GB" "${dest}"
  for domain in bac120 ar53; do
    tsv="${dest}/${domain}_metadata_r${GTDB_RELEASE}.tsv"
    if [[ -s "${tsv}" ]]; then
      log_info "present: $(basename "${tsv}")"
      continue
    fi
    rel="${domain}_metadata_r${GTDB_RELEASE}.tsv.gz"
    fetch "${rel_base}/${rel}" "${tsv}.gz"
    expected="$(gtdb_md5 "${rel}")"
    if [[ -n "${expected}" ]]; then verify_md5 "${tsv}.gz" "${expected}"; else log_warn "no published checksum for ${rel}"; fi
    gunzip -f "${tsv}.gz"
  done

  announce "NCBI taxdump (for genus..phylum TaxId resolution)" "~0.5 GB" "${dest}/ncbi_taxdump"
  if [[ -s "${dest}/ncbi_taxdump/names.dmp" ]]; then
    log_info "present, skipping"
  else
    td="${dest}/downloads/taxdump.tar.gz"
    fetch "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz" "${td}"
    fetch "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5" "${td}.md5"
    verify_md5 "${td}" "$(awk '{ print $1 }' "${td}.md5")"
    mkdir -p "${dest}/ncbi_taxdump"
    tar -xzf "${td}" -C "${dest}/ncbi_taxdump" names.dmp nodes.dmp
    rm -f "${td}" "${td}.md5"
  fi
  # The Mash sketch is built by GTDB-Tk on first use and cached in mash_db/.
fi

# ----------------------------------------------------------------------------
# DRAM
# ----------------------------------------------------------------------------
if wanted dram; then
  dest="${DB_ROOT}/dram"
  announce "DRAM databases" "~500 GB" "${dest}"
  DRAM_OPTS=(--output_dir "${dest}" --threads "${THREADS}")
  if is_true "${SKIP_UNIREF}"; then
    # UniRef90 is the long pole (download + MMseqs2 indexing). Without it,
    # genes lacking KEGG/Pfam/dbCAN/MEROPS hits get no UniRef annotation and
    # DRAM's distill has less to work with when KEGG is also absent.
    log_warn "--skip-uniref: UniRef90 will not be available to DRAM annotate"
    DRAM_OPTS+=(--skip_uniref)
  fi
  if [[ -n "${DRAM_KEGG_LOC:-}" ]]; then
    require_file "${DRAM_KEGG_LOC}"
    DRAM_OPTS+=(--kegg_loc "${DRAM_KEGG_LOC}")
  else
    log_info "DRAM_KEGG_LOC not set; KEGG (licensed, optional) is skipped"
  fi
  if [[ -d "${dest}" ]] && [[ -n "$(ls -A "${dest}" 2>/dev/null)" ]]; then
    log_warn "${dest} is not empty; DRAM-setup.py will reuse/overwrite its contents"
  fi
  # DRAM verifies downloads internally and records the paths in its env config.
  run_in_env dram DRAM-setup.py prepare_databases "${DRAM_OPTS[@]}"
  run_in_env dram DRAM-setup.py print_config
fi

# ----------------------------------------------------------------------------
# MetaPhlAn
# ----------------------------------------------------------------------------
if wanted metaphlan; then
  dest="${DB_ROOT}/metaphlan"
  announce "MetaPhlAn database" "~20 GB" "${dest}"
  mkdir -p "${dest}"
  # metaphlan --install checks the md5 of the downloaded bundle itself.
  run_in_env metaphlan metaphlan --install --db_dir "${dest}" --nproc "${THREADS}"
fi

# ----------------------------------------------------------------------------
# mettannotator
# ----------------------------------------------------------------------------
# mettannotator has no standalone download script: it fills --dbs with
# whatever is missing on first run. To pre-stage, run it once on a single
# small genome (Mycoplasma genitalium G37, NCBI TaxId 2097) so the first real
# run is not also a 182 GB download run. Containers are cached too.
if wanted mettannotator; then
  dest="${METTANNOTATOR_DB}"
  announce "mettannotator databases + containers" "~182 GB" "${dest}"
  load_env_module "${NEXTFLOW_MODULE}"
  load_env_module "${APPTAINER_MODULE}"
  require_cmd nextflow
  warm="${DB_ROOT}/mettannotator_warmup"
  mkdir -p "${dest}" "${NXF_SINGULARITY_CACHEDIR}" "${warm}"
  export NXF_SINGULARITY_CACHEDIR
  export NXF_WORK="${warm}/work"
  export NXF_HOME="${NXF_HOME:-${DB_ROOT}/nextflow_home}"
  fetch "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000908.2&rettype=fasta&retmode=text" \
    "${warm}/MGEN_G37.fa"
  require_file "${warm}/MGEN_G37.fa"
  printf 'prefix,assembly,taxid\nMGEN_G37,%s,2097\n' "${warm}/MGEN_G37.fa" > "${warm}/samplesheet.csv"
  (
    cd "${warm}"
    nextflow run ebi-metagenomics/mettannotator \
      -r "${METTANNOTATOR_VERSION}" \
      -profile singularity \
      --input "${warm}/samplesheet.csv" \
      --outdir "${warm}/results" \
      --dbs "${dest}"
  )
  log_info "mettannotator databases staged in ${dest}; ${warm} can be deleted"
fi

summary_block "download_databases: done" \
  "DB_ROOT" "${DB_ROOT}" \
  "Requested" "${ONLY:-all}" \
  "DRAM UniRef90 skipped" "${SKIP_UNIREF}"
