#!/usr/bin/env bash
#SBATCH --job-name=mag11_metapop_prep
#SBATCH --output=mag11_metapop_prep_%j.out
#SBATCH --error=mag11_metapop_prep_%j.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 11_metapop_prep.sh
# Build the shared reference that every sample is remapped to for MetaPop.
#
# Contigs are reheaded "<MAG>_<contig>". The reheaded per-genome FASTAs are
# written to metapop/reference/ and that directory (not the dRep directory)
# is MetaPop's --reference: MetaPop matches BAM contig names against the
# reference FASTAs, so the names must be identical in both.
#
# Outputs:
#   metapop/reference/<MAG>.fa            reheaded per-genome FASTAs
#   metapop/combined_mags.fasta           concatenation of the above
#   metapop/idx/combined_mags.*.bt2[l]    bowtie2 index
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "11_metapop_prep"
exit_if_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"

REPS_DIR="$(drep_reps_dir)"
METAPOP_DIR="${OUT_DIR}/metapop"
REF_DIR="${METAPOP_DIR}/reference"
COMBINED="${METAPOP_DIR}/combined_mags.fasta"
IDX="${METAPOP_DIR}/idx/combined_mags"

require_dir "${REPS_DIR}"
N_REPS="$(count_files "${REPS_DIR}" '*.fa')"
[[ "${N_REPS}" -gt 0 ]] || die "no representative genomes in ${REPS_DIR}"

# ----------------------------------------------------------------------------
# Rehead and concatenate
# ----------------------------------------------------------------------------
rm -rf "${REF_DIR}" "${METAPOP_DIR}/idx"
mkdir -p "${REF_DIR}" "${METAPOP_DIR}/idx"
: > "${COMBINED}"
while IFS= read -r fa; do
  mag="$(basename "${fa}" .fa)"
  awk -v m="${mag}" '/^>/ { print ">" m "_" substr($1, 2); next } { print }' \
    "${fa}" > "${REF_DIR}/${mag}.fa"
  cat "${REF_DIR}/${mag}.fa" >> "${COMBINED}"
done < <(find "${REPS_DIR}" -maxdepth 1 -type f -name '*.fa' | sort)

N_CONTIGS="$(awk '/^>/ { n++ } END { print n + 0 }' "${COMBINED}")"
N_UNIQUE="$(grep '^>' "${COMBINED}" | sort -u | wc -l | tr -d ' ')"
[[ "${N_CONTIGS}" -eq "${N_UNIQUE}" ]] || die "combined reference has duplicate contig names (${N_CONTIGS} contigs, ${N_UNIQUE} unique)"
log_info "combined reference: ${N_REPS} genomes, ${N_CONTIGS} contigs"

# ----------------------------------------------------------------------------
# bowtie2 index
# ----------------------------------------------------------------------------
log_info "bowtie2-build"
run_in_env metapop bowtie2-build --threads "${THREADS}" "${COMBINED}" "${IDX}" > /dev/null

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
summary_block "11_metapop_prep: summary" \
  "Genomes in reference" "${N_REPS}" \
  "Total contigs" "${N_CONTIGS}" \
  "Reference directory" "${REF_DIR}" \
  "Combined FASTA" "${COMBINED}" \
  "bowtie2 index" "${IDX}"

mark_done "${MODULE_NAME}" "${ALL_SAMPLES_KEY}"
