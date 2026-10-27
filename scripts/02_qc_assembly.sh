#!/usr/bin/env bash
#SBATCH --job-name=mag02_qc_assembly
#SBATCH --output=mag02_qc_assembly_%A_%a.out
#SBATCH --error=mag02_qc_assembly_%A_%a.err
set -Eeuo pipefail

# ----------------------------------------------------------------------------
# 02_qc_assembly.sh  (ARRAY: one task per line of SAMPLE_LIST)
# fastp -> host removal -> MicrobeCensus -> MetaPhlAn -> fastp dedup ->
# metaSPAdes -> QUAST
#
# fastp runs twice, producing two read sets from one host-removal pass:
#   QC reads     trimmed + host-removed, NOT deduplicated. Used for profiling
#                (MicrobeCensus, MetaPhlAn) and everything coverage-based
#                downstream (03 binning depths, 12 MetaPop remapping).
#   dedup reads  the QC reads after a second, dedup-only fastp pass on the
#                lane-merged pair. Used only for assembly.
#
# Outputs:
#   fastp/<S>_L<i>.{json,html}                  QC pass, per library
#   fastp/<S>.dedup.{json,html}                 dedup pass
#   qcreads/<S>_L<i>_R{1,2}.fastq.gz            (when KEEP_QC_READS=true)
#   qcreads/dedup/<S>_R{1,2}.fastq.gz           (when KEEP_QC_READS=true)
#   microbecensus/<S>.microbecensus
#   metaphlan/<S>.metaphlan
#   assemblies/<S>_scaffolds.fasta
#   assemblies/assembly_stats/<S>/               (QUAST)
# ----------------------------------------------------------------------------

if [[ -z "${MAG_SLURM_ROOT:-}" ]]; then
  MAG_SLURM_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi
# shellcheck source-path=SCRIPTDIR source=../lib/common.sh
source "${MAG_SLURM_ROOT}/lib/common.sh"
init_module "02_qc_assembly"

SAMPLE="$(sample_for_task "${SLURM_ARRAY_TASK_ID:-${1:-}}")"
exit_if_done "${MODULE_NAME}" "${SAMPLE}"
start_sample_log "${SAMPLE}"

LIB_COUNT="$(library_count "${SAMPLE}")"
MEM_GB="$(mem_gb "${MEM_QC_ASSEMBLY}")"
setup_workdir

# ----------------------------------------------------------------------------
# 1. Stage reads
# ----------------------------------------------------------------------------
# R1[i] / R2[i] always hold the current read files for library i; each step
# below replaces them with its own output.
mkdir -p "${WORKDIR}/reads"
R1=()
R2=()
RAW_FILES=""
for (( i = 1; i <= LIB_COUNT; i++ )); do
  col_f=$(( 2 * i ))
  col_r=$(( 2 * i + 1 ))
  src_f="${READ_DIR}/$(manifest_field "${SAMPLE}" "${col_f}")"
  src_r="${READ_DIR}/$(manifest_field "${SAMPLE}" "${col_r}")"
  require_file "${src_f}" "${src_r}"
  R1[i]="${WORKDIR}/reads/${SAMPLE}_L${i}_R1.fastq.gz"
  R2[i]="${WORKDIR}/reads/${SAMPLE}_L${i}_R2.fastq.gz"
  cp "${src_f}" "${R1[i]}"
  cp "${src_r}" "${R2[i]}"
  RAW_FILES="${RAW_FILES:+${RAW_FILES},}${R1[i]},${R2[i]}"
done
log_info "staged ${LIB_COUNT} librar$( [[ "${LIB_COUNT}" -eq 1 ]] && echo y || echo ies)"

# ----------------------------------------------------------------------------
# 2. fastp, pass 1: trimming and filtering (no dedup)
# ----------------------------------------------------------------------------
if is_true "${RUN_FASTP}"; then
  log_info "fastp pass 1 (QC)"
  mkdir -p "${WORKDIR}/fastp" "${OUT_DIR}/fastp"
  FASTP_OPTS=(--detect_adapter_for_pe --cut_front --cut_tail
    --cut_window_size "${FASTP_CUT_WINDOW_SIZE}"
    --cut_mean_quality "${FASTP_CUT_MEAN_QUALITY}"
    --length_required "${FASTP_LENGTH_REQUIRED}")
  if is_true "${FASTP_TRIM_POLY_G}"; then
    FASTP_OPTS+=(--trim_poly_g)
  fi
  # No --dedup here: duplicate removal would distort the coverage that
  # MetaBAT2 and MetaPop use (legacy defect #15). Dedup is pass 2, for
  # assembly only.
  for (( i = 1; i <= LIB_COUNT; i++ )); do
    out1="${WORKDIR}/fastp/${SAMPLE}_L${i}_R1.fastq.gz"
    out2="${WORKDIR}/fastp/${SAMPLE}_L${i}_R2.fastq.gz"
    run_in_env fastp fastp \
      --in1 "${R1[i]}" --in2 "${R2[i]}" \
      --out1 "${out1}" --out2 "${out2}" \
      "${FASTP_OPTS[@]}" \
      --json "${OUT_DIR}/fastp/${SAMPLE}_L${i}.json" \
      --html "${OUT_DIR}/fastp/${SAMPLE}_L${i}.html" \
      --thread "${THREADS}"
    R1[i]="${out1}"
    R2[i]="${out2}"
  done
fi

# ----------------------------------------------------------------------------
# 3. Host removal
# ----------------------------------------------------------------------------
if is_true "${RUN_HOST_REMOVAL}" && [[ "${HOST_SPECIES}" != "none" ]]; then
  HOST_INDEX="$(resolve_host_index)"
  log_info "host removal against ${HOST_INDEX}"
  mkdir -p "${WORKDIR}/host"
  for (( i = 1; i <= LIB_COUNT; i++ )); do
    # Alignments are discarded (-S /dev/null); only the concordantly unaligned
    # pairs are kept. The legacy code wrote and then deleted a full SAM
    # (legacy defect #10).
    run_in_env fastp bowtie2 \
      --very-sensitive --dovetail \
      -p "${THREADS}" \
      -x "${HOST_INDEX}" \
      -1 "${R1[i]}" -2 "${R2[i]}" \
      --un-conc-gz "${WORKDIR}/host/${SAMPLE}_L${i}_R%.fastq.gz" \
      -S /dev/null
    R1[i]="${WORKDIR}/host/${SAMPLE}_L${i}_R1.fastq.gz"
    R2[i]="${WORKDIR}/host/${SAMPLE}_L${i}_R2.fastq.gz"
    require_file "${R1[i]}" "${R2[i]}"
  done
else
  log_info "host removal skipped (RUN_HOST_REMOVAL=${RUN_HOST_REMOVAL}, HOST_SPECIES=${HOST_SPECIES})"
fi

# ----------------------------------------------------------------------------
# 4. Keep post-QC reads
# ----------------------------------------------------------------------------
# 03_binning.sh and 12_metapop_remap.sh depend on this exact naming.
if is_true "${KEEP_QC_READS}"; then
  mkdir -p "${QC_READ_DIR}"
  for (( i = 1; i <= LIB_COUNT; i++ )); do
    cp "${R1[i]}" "${QC_READ_DIR}/${SAMPLE}_L${i}_R1.fastq.gz"
    cp "${R2[i]}" "${QC_READ_DIR}/${SAMPLE}_L${i}_R2.fastq.gz"
  done
  log_info "post-QC reads copied to ${QC_READ_DIR}"
fi

# ----------------------------------------------------------------------------
# 5. Merge libraries
# ----------------------------------------------------------------------------
# Concatenated gzip members form a valid gzip stream, so `cat` is sufficient.
MERGED_R1="${WORKDIR}/reads/${SAMPLE}_R1.fastq.gz"
MERGED_R2="${WORKDIR}/reads/${SAMPLE}_R2.fastq.gz"
if [[ "${LIB_COUNT}" -eq 2 ]]; then
  log_info "merging 2 libraries"
  cat "${R1[1]}" "${R1[2]}" > "${MERGED_R1}"
  cat "${R2[1]}" "${R2[2]}" > "${MERGED_R2}"
else
  mv "${R1[1]}" "${MERGED_R1}"
  mv "${R2[1]}" "${MERGED_R2}"
fi
MERGED_SINGLE="${WORKDIR}/reads/${SAMPLE}_merged.fastq.gz"
if is_true "${RUN_METAPHLAN}"; then
  cat "${MERGED_R1}" "${MERGED_R2}" > "${MERGED_SINGLE}"
fi

# ----------------------------------------------------------------------------
# 6. MicrobeCensus
# ----------------------------------------------------------------------------
if is_true "${RUN_MICROBECENSUS}"; then
  # Average genome size estimation is read-length sensitive; "raw" (the
  # default) runs on untrimmed reads (legacy defect #16).
  if [[ "${MICROBECENSUS_INPUT}" == "raw" ]]; then
    MC_INPUT="${RAW_FILES}"
  else
    MC_INPUT="${MERGED_R1},${MERGED_R2}"
  fi
  log_info "MicrobeCensus on ${MICROBECENSUS_INPUT} reads"
  mkdir -p "${OUT_DIR}/microbecensus"
  run_in_env microbecensus run_microbe_census.py \
    -v -t "${THREADS}" \
    "${MC_INPUT}" \
    "${OUT_DIR}/microbecensus/${SAMPLE}.microbecensus"
fi

# ----------------------------------------------------------------------------
# 7. MetaPhlAn
# ----------------------------------------------------------------------------
if is_true "${RUN_METAPHLAN}"; then
  log_info "MetaPhlAn"
  mkdir -p "${OUT_DIR}/metaphlan" "${WORKDIR}/metaphlan_tmp"
  run_in_env metaphlan metaphlan \
    "${MERGED_SINGLE}" \
    --input_type fastq \
    -t rel_ab_w_read_stats \
    --db_dir "${DB_ROOT}/metaphlan" \
    --sample_id "${SAMPLE}" \
    --nproc "${THREADS}" \
    --mapout "${WORKDIR}/${SAMPLE}.mapout.bz2" \
    --tmp_dir "${WORKDIR}/metaphlan_tmp" \
    --output_file "${OUT_DIR}/metaphlan/${SAMPLE}.metaphlan"
  rm -f "${MERGED_SINGLE}"
fi

# ----------------------------------------------------------------------------
# 8. fastp, pass 2: deduplication for assembly
# ----------------------------------------------------------------------------
# Runs on the lane-merged QC pair so duplicates shared across lanes of the same
# library are also removed. Every other fastp operation is disabled (-A
# adapters, -G poly-G, -Q quality filter, -L length filter): pass 1 already
# did them, and this pass must only drop duplicates.
ASM_R1="${MERGED_R1}"
ASM_R2="${MERGED_R2}"
if is_true "${RUN_ASSEMBLY}" && is_true "${FASTP_DEDUP_ASSEMBLY}"; then
  log_info "fastp pass 2 (dedup for assembly)"
  mkdir -p "${WORKDIR}/dedup" "${OUT_DIR}/fastp"
  ASM_R1="${WORKDIR}/dedup/${SAMPLE}_R1.fastq.gz"
  ASM_R2="${WORKDIR}/dedup/${SAMPLE}_R2.fastq.gz"
  run_in_env fastp fastp \
    --in1 "${MERGED_R1}" --in2 "${MERGED_R2}" \
    --out1 "${ASM_R1}" --out2 "${ASM_R2}" \
    --dedup \
    -A -G -Q -L \
    --json "${OUT_DIR}/fastp/${SAMPLE}.dedup.json" \
    --html "${OUT_DIR}/fastp/${SAMPLE}.dedup.html" \
    --thread "${THREADS}"
  require_file "${ASM_R1}" "${ASM_R2}"
  if is_true "${KEEP_QC_READS}"; then
    mkdir -p "${QC_READ_DIR}/dedup"
    cp "${ASM_R1}" "${QC_READ_DIR}/dedup/${SAMPLE}_R1.fastq.gz"
    cp "${ASM_R2}" "${QC_READ_DIR}/dedup/${SAMPLE}_R2.fastq.gz"
    log_info "deduplicated reads copied to ${QC_READ_DIR}/dedup"
  fi
elif is_true "${RUN_ASSEMBLY}"; then
  log_info "FASTP_DEDUP_ASSEMBLY=false: assembling the non-deduplicated QC reads"
fi

# ----------------------------------------------------------------------------
# 9. metaSPAdes
# ----------------------------------------------------------------------------
ASSEMBLY="${OUT_DIR}/assemblies/${SAMPLE}_scaffolds.fasta"
if is_true "${RUN_ASSEMBLY}"; then
  log_info "metaSPAdes on $( [[ "${ASM_R1}" == "${MERGED_R1}" ]] && echo "QC" || echo "deduplicated") reads (${THREADS} threads, ${MEM_GB} GB)"
  mkdir -p "${WORKDIR}/spades_tmp" "${OUT_DIR}/assemblies"
  # metaSPAdes accepts exactly one paired-end library ("--pe-1 1 ... --pe-1 2"
  # is rejected in --meta mode). Two lanes of one sample are the same library,
  # so a single lane-merged pair is passed in both cases.
  run_in_env spades metaspades.py \
    --threads "${THREADS}" \
    --memory "${MEM_GB}" \
    --tmp-dir "${WORKDIR}/spades_tmp" \
    -1 "${ASM_R1}" -2 "${ASM_R2}" \
    -o "${WORKDIR}/spades"
  require_file "${WORKDIR}/spades/scaffolds.fasta"
  cp "${WORKDIR}/spades/scaffolds.fasta" "${ASSEMBLY}"
  cp "${WORKDIR}/spades/spades.log" "${OUT_DIR}/assemblies/${SAMPLE}_spades.log" 2>/dev/null || true
fi

# ----------------------------------------------------------------------------
# 10. QUAST
# ----------------------------------------------------------------------------
if is_true "${RUN_QUAST}"; then
  require_file "${ASSEMBLY}"
  log_info "QUAST"
  mkdir -p "${OUT_DIR}/assemblies/assembly_stats/${SAMPLE}"
  run_in_env quast metaquast.py \
    --threads "${THREADS}" \
    --max-ref-number 0 \
    --output-dir "${OUT_DIR}/assemblies/assembly_stats/${SAMPLE}" \
    "${ASSEMBLY}"
fi

# ----------------------------------------------------------------------------
# Summary
# ----------------------------------------------------------------------------
N_CONTIGS="NA"
if [[ -s "${ASSEMBLY}" ]]; then
  N_CONTIGS="$(awk '/^>/ { n++ } END { print n + 0 }' "${ASSEMBLY}")"
fi
mark_done "${MODULE_NAME}" "${SAMPLE}"
summary_block "02_qc_assembly: ${SAMPLE}" \
  "Libraries" "${LIB_COUNT}" \
  "fastp / host removal" "${RUN_FASTP} / ${RUN_HOST_REMOVAL} (${HOST_SPECIES})" \
  "QC reads kept" "${KEEP_QC_READS}" \
  "Assembly reads deduplicated" "$( is_true "${RUN_ASSEMBLY}" && echo "${FASTP_DEDUP_ASSEMBLY}" || echo "n/a")" \
  "MicrobeCensus / MetaPhlAn" "${RUN_MICROBECENSUS} / ${RUN_METAPHLAN}" \
  "Assembly" "${ASSEMBLY}" \
  "Scaffolds" "${N_CONTIGS}"
