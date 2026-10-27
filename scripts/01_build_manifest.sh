#!/usr/bin/env bash
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 01_build_manifest.sh
# Interactive helper (NOT a SLURM job). Scans READ_DIR for paired FASTQs and
# writes MANIFEST (with header) and SAMPLE_LIST (headerless).
#
# Usage:
#   scripts/01_build_manifest.sh [--read-dir DIR] [--pattern PATTERN] [--force]
#
# PATTERN is a filename template with three placeholders:
#   *        the SampleID
#   {lane}   the lane / library number (digits); omit for single-lane data
#   {1,2}    the mate (1 or 2)
# Default: '*_L{lane}_R{1,2}.fastq.gz'
# Example for SRA fasterq-dump output: '*_{1,2}.fastq.gz'
# Only one '*' is allowed; Illumina '*_S*_L001_R1_001' names must be renamed.
# ----------------------------------------------------------------------------

MAG_SLURM_ROOT="${MAG_SLURM_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"

MODULE_NAME="01_build_manifest"
install_traps

usage() {
  sed -n '5,19p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

# ----------------------------------------------------------------------------
# Arguments
# ----------------------------------------------------------------------------
PATTERN='*_L{lane}_R{1,2}.fastq.gz'
READ_DIR_ARG=""
FORCE=false
while [[ $# -gt 0 ]]; do
  case "$1" in
    --read-dir) READ_DIR_ARG="${2:?--read-dir needs a value}"; shift 2 ;;
    --pattern) PATTERN="${2:?--pattern needs a value}"; shift 2 ;;
    --config) MAG_SLURM_CONFIG="${2:?--config needs a value}"; export MAG_SLURM_CONFIG; shift 2 ;;
    --force) FORCE=true; shift ;;
    -h | --help) usage; exit 0 ;;
    *) die "unknown argument: $1 (see --help)" ;;
  esac
done

load_config
validate_config --skip-manifest
if [[ -n "${READ_DIR_ARG}" ]]; then
  READ_DIR="${READ_DIR_ARG}"
fi
require_dir "${READ_DIR}"

if [[ -e "${MANIFEST}" || -e "${SAMPLE_LIST}" ]] && ! is_true "${FORCE}"; then
  die "${MANIFEST} or ${SAMPLE_LIST} already exists; pass --force to overwrite"
fi

# ----------------------------------------------------------------------------
# Translate the pattern into an anchored regex
# ----------------------------------------------------------------------------
[[ "${PATTERN}" == *'{1,2}'* ]] || die "pattern must contain the mate placeholder {1,2}: ${PATTERN}"
n_stars="$(printf '%s' "${PATTERN}" | tr -cd '*' | wc -c | tr -d ' ')"
[[ "${n_stars}" -eq 1 ]] || die "pattern must contain exactly one '*' (the SampleID): ${PATTERN}"

# Escape regex metacharacters, then substitute placeholders. The order in which
# placeholders appear determines the capture group numbers.
# shellcheck disable=SC2016  # literal regex metacharacters, not expansions
regex="$(printf '%s' "${PATTERN}" \
  | sed -e 's/[].[^$()+?|\\]/\\&/g' \
        -e 's/{lane}/@LANE@/' -e 's/{1,2}/@MATE@/' -e 's/\*/@SAMPLE@/')"
regex="^${regex}$"
group_order="$(printf '%s' "${regex}" | grep -o '@[A-Z]*@' | tr -d '@' | tr '\n' ' ')"
regex="${regex//@SAMPLE@/(.+)}"
regex="${regex//@LANE@/([0-9]+)}"
regex="${regex//@MATE@/([12])}"
log_info "pattern ${PATTERN} -> regex ${regex} (groups: ${group_order})"

# ----------------------------------------------------------------------------
# Scan READ_DIR
# ----------------------------------------------------------------------------
TMP_LIST="$(mktemp)"
UNMATCHED="$(mktemp)"
CLEANUP_PATHS="${TMP_LIST} ${UNMATCHED}"

shopt -s nullglob
for path in "${READ_DIR}"/*; do
  [[ -f "${path}" ]] || continue
  name="$(basename "${path}")"
  if [[ "${name}" =~ ${regex} ]]; then
    sample=""; lane="1"; mate=""
    g=1
    for kind in ${group_order}; do
      case "${kind}" in
        SAMPLE) sample="${BASH_REMATCH[${g}]}" ;;
        LANE) lane="${BASH_REMATCH[${g}]}" ;;
        MATE) mate="${BASH_REMATCH[${g}]}" ;;
      esac
      g=$(( g + 1 ))
    done
    printf '%s\t%s\t%s\t%s\n' "${sample}" "${lane}" "${mate}" "${name}" >> "${TMP_LIST}"
  elif [[ "${name}" == *.fastq.gz || "${name}" == *.fq.gz || "${name}" == *.fastq || "${name}" == *.fq ]]; then
    printf '%s\n' "${name}" >> "${UNMATCHED}"
  fi
done
shopt -u nullglob

[[ -s "${TMP_LIST}" ]] || die "no files in ${READ_DIR} match ${PATTERN}"
if [[ -s "${UNMATCHED}" ]]; then
  log_warn "$(count_lines "${UNMATCHED}") FASTQ file(s) in ${READ_DIR} do not match the pattern and are NOT in the manifest:"
  head -n 20 "${UNMATCHED}" | sed 's/^/    /' >&2
fi

# ----------------------------------------------------------------------------
# Group into samples and validate pairing
# ----------------------------------------------------------------------------
# Unpaired reads are an error, not a silent drop. Samples with more than two
# lanes cannot be represented in the 5-column manifest and are also an error.
MANIFEST_TMP="${MANIFEST}.tmp"
if ! sort -t "$(printf '\t')" -k1,1 -k2,2n -k3,3n "${TMP_LIST}" | awk -F'\t' '
  {
    key = $1 SUBSEP $2
    if ((key SUBSEP $3) in seen) { printf "ERROR: duplicate file for sample %s lane %s mate %s: %s\n", $1, $2, $3, $4 > "/dev/stderr"; bad = 1 }
    seen[key SUBSEP $3] = $4
    if (!($1 in nlanes)) { order[++ns] = $1; nlanes[$1] = 0 }
    if (!(key in lane_seen)) { lane_seen[key] = 1; nlanes[$1]++; lanes[$1, nlanes[$1]] = $2 }
  }
  END {
    print "SampleID\tLib1_F\tLib1_R\tLib2_F\tLib2_R"
    for (i = 1; i <= ns; i++) {
      s = order[i]
      if (nlanes[s] > 2) { printf "ERROR: sample %s has %d lanes; the manifest supports at most 2 (concatenate lanes first)\n", s, nlanes[s] > "/dev/stderr"; bad = 1; continue }
      row = s
      for (j = 1; j <= nlanes[s]; j++) {
        l = lanes[s, j]
        f1 = ((s SUBSEP l SUBSEP 1) in seen) ? seen[s SUBSEP l SUBSEP 1] : ""
        f2 = ((s SUBSEP l SUBSEP 2) in seen) ? seen[s SUBSEP l SUBSEP 2] : ""
        if (f1 == "" || f2 == "") {
          printf "ERROR: unpaired read for sample %s lane %s: R1=%s R2=%s\n", s, l, (f1 == "" ? "MISSING" : f1), (f2 == "" ? "MISSING" : f2) > "/dev/stderr"
          bad = 1
        }
        row = row "\t" f1 "\t" f2
      }
      if (nlanes[s] == 1) row = row "\t\t"
      print row
    }
    exit bad
  }' > "${MANIFEST_TMP}"; then
  rm -f "${MANIFEST_TMP}"
  die "manifest validation failed (see errors above); nothing was written"
fi

# ----------------------------------------------------------------------------
# SampleID sanity
# ----------------------------------------------------------------------------
n_bad_ids=0
while IFS= read -r sid; do
  if [[ "${sid}" =~ [[:space:]] || "${sid}" == */* || "${sid}" == -* ]]; then
    log_warn "SampleID '${sid}' contains whitespace, '/', or a leading '-'; this breaks downstream tools"
    n_bad_ids=$(( n_bad_ids + 1 ))
  fi
done < <(tail -n +2 "${MANIFEST_TMP}" | cut -f 1)

mv "${MANIFEST_TMP}" "${MANIFEST}"
# samples.txt is headerless by definition (legacy defect #6).
tail -n +2 "${MANIFEST}" | cut -f 1 > "${SAMPLE_LIST}"

n_samples="$(count_lines "${SAMPLE_LIST}")"
n_two_lane="$(awk -F'\t' 'NR > 1 && $4 != "" { n++ } END { print n + 0 }' "${MANIFEST}")"
[[ "${n_samples}" -eq $(( $(count_lines "${MANIFEST}") - 1 )) ]] \
  || die "internal error: SAMPLE_LIST and MANIFEST row counts differ"

summary_block "Manifest written" \
  "Read directory" "${READ_DIR}" \
  "Pattern" "${PATTERN}" \
  "Samples" "${n_samples}" \
  "Samples with 2 libraries" "${n_two_lane}" \
  "Unmatched FASTQs (ignored)" "$( [[ -s "${UNMATCHED}" ]] && count_lines "${UNMATCHED}" || echo 0)" \
  "SampleIDs with bad characters" "${n_bad_ids}" \
  "Manifest" "${MANIFEST}" \
  "Sample list" "${SAMPLE_LIST}"
