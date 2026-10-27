# Main Pipeline

> **SUPERSEDED — DO NOT EXECUTE.**
> This is the prior hand-written pipeline, retained verbatim for provenance and as the
> reference point for the bug-fix checklist in the build spec. It contains known defects
> (see §9 of the build spec) and hardcoded paths specific to one user and one cluster.
> It is excluded from linting and CI. Destination: `docs/legacy/original_pipeline.md`.

---

Reads were downloaded from the SRA with `fasterq-dump` under the BioProject PRJNA662107.

## SRA Download

```{bash, eval = F}
#!/bin/bash
#SBATCH --job-name=sra_fasterq
#SBATCH --output=logs/sra_fasterq_%A_%a.out
#SBATCH --error=logs/sra_fasterq_%A_%a.err
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --array=1-523%10
#SBATCH --partition=sla-prio
#SBATCH --account=one_sc_default
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=stc5395@psu.edu

# -----------------------------
# User settings
# -----------------------------

ACCESSION_LIST="/storage/home/stc5395/scratch/BordensteinVegVirome/SRR_Acc_List.txt"
OUTDIR="/storage/home/stc5395/scratch/BordensteinVegVirome/raw_reads"
TMPDIR="/storage/home/stc5395/scratch/BordensteinVegVirome/sra_tmp/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

THREADS="${SLURM_CPUS_PER_TASK}"

source /storage/home/stc5395/work/miniconda3/etc/profile.d/conda.sh

# -----------------------------
# Setup
# -----------------------------

mkdir -p "$OUTDIR" "$TMPDIR" logs

ACCESSION=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$ACCESSION_LIST")

if [[ -z "${ACCESSION}" ]]; then
    echo "No accession found for array task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Downloading accession: ${ACCESSION}"
echo "Output directory: ${OUTDIR}"
echo "Temp directory: ${TMPDIR}"

# -----------------------------
# Download and convert
# -----------------------------

conda run -n sratool_env fasterq-dump "$ACCESSION" \
    --outdir "$OUTDIR" \
    --temp "$TMPDIR" \
    --threads "$THREADS" \
    --split-files \
    --progress

# Optional compression
gzip -f "$OUTDIR/${ACCESSION}"*.fastq

# Cleanup
rm -rf "$TMPDIR"

echo "Finished ${ACCESSION}"
```

Next reformat the file names such that they are compatible with the lab metagenomic assembly script.

```{bash, eval = F}
cd /storage/home/stc5395/scratch/BordensteinVegVirome/raw_reads

for f in SRR*_1.fastq.gz SRR*_2.fastq.gz; do
  new=$(echo "$f" | sed -E 's/_1\.fastq\.gz/_L1_R1.fastq.gz/; s/_2\.fastq\.gz/_L1_R2.fastq.gz/')
  mv "$f" "$new"
done
```

Generate the manifest file with the following headers: SampleID	Lib1_F Lib1_R

```{r, eval = F}
read_csv("/home/stc5395/projects/BordensteinVegVirome/data/runinfo.csv") %>% 
  dplyr::select(Run) %>%
  dplyr::rename(SampleID = Run) %>%
  mutate(
    Lib1_F = paste0(SampleID, "_L1_R1.fastq.gz"),
    Lib1_R = paste0(SampleID, "_L1_R2.fastq.gz")
  ) %>%
  write_tsv("/home/stc5395/projects/BordensteinVegVirome/data/VegViromeAssemblyManifest.tsv")

read_csv("/home/stc5395/projects/BordensteinVegVirome/data/runinfo.csv") %>% 
  dplyr::select(Run) %>%
  write_tsv("/home/stc5395/projects/BordensteinVegVirome/data/VegViromeAssemblySample.tsv")
```

## Fastp, MetaSPAdes, MetaPHlAn, and QUAST

```{bash, eval = F}
#!/usr/bin/env bash
#SBATCH --job-name VegVirAssembly
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=256GB
#SBATCH --array=1-523 # the number of input files for parallel processing. This would do 130 samples from the manifest starting with the first. Note that the headerline is not counted.
#SBATCH --output assembly_pipeline.log
#SBATCH --error assembly_pipeline.err
#SBATCH --partition=sla-prio
#SBATCH --account=one_sc_default
#SBATCH --time=48:00:00 #change as needed
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user=stc5395@psu.edu

######################################################################################################################################################################################
# Metagenome Pipeline v0.011 JB
# Notes:
# -if interactive session is needed: salloc --nodes=1 --ntasks=4 --mem=32G --time=03:00:00
# --array=1-2 which lines of the file to go through. This would do the first two samples (note the headerline is not counted)
# -on ROAR large files should go on the scratch directory, but will be cleared off after 30 days. It would be ideal to stage the results here but return the final files to work or home
# -Can also use account=open amd remove partition line, as SBATCH --account=open
# - this script can also assemble; however, may required longer than 48h. Recommend to set RunMicrobeCensus, RunMetaphlan, and RunHumann to false for assembly.
######################################################################################################################################################################################

######################################################################################################################################################################################
# Change log:
# - 30 Apr 2025 - initial setup and debugging
# - 28 July 2025 - changed handling of postqc reads, added explicit polyG removal assuming all data will be novaseq/nextseq, although this would have been occurring by default
# - 24 June 2026 -  Changed location to reflect personal directory structure, utilized individualized conda environments for each command to avoid conflicts. Added quast for assembly qc
#
#
#
#
######################################################################################################################################################################################

#/storage/work/jpb6325/conda_roar is install location
#created humann_30Apr2025 environment with humann, kneaddata, fastp, and metaspades. Microbecensus and MetaSpades  are in own environments due to python conflict

# Load humann environment from JB's work directory. This is humann 3.9
#humann_config --update database_folders nucleotide /storage/work/jpb6325/humann_databases/chocophlan/
#humann_config --update database_folders protein /storage/work/jpb6325/humann_databases/uniref/
#humann_config --update database_folders utility_mapping /storage/work/jpb6325/humann_databases/utility_mapping/

######################################################################################################################################################################################
# Tasks -> set true or false as desired
RunFastP=true # do QC on reads
RemoveHost=true #Note: adjust as needed below, defaulting to human host
HostSpecies="human" # can be mouse or human which controls reference genome for removing host reads
KeepQCReads=true # should the QC reads be kept? This might not be preferable for storage reasons, depends on application. Be sure to store the reads somewhere with sufficient storage (i.e. not home directory)
RunMicrobeCensus=true # calculate genome equivalent calculations
RunMetaphlan=true # Should metaphlan be run? Note: do not have this be true if humann is also running as they will be redundant 
RunHumann=false # Should HUMANN be run?
RunMetaSpades=true # Should assembly be done?

######################################################################################################################################################################################
#Input/Putput parameters -> be sure to use absolute file paths for cluster
OutDir="/storage/home/stc5395/scratch/BordensteinVegVirome"
Reads="/storage/home/stc5395/scratch/BordensteinVegVirome/raw_reads" # a directory containing all fastq files to be processed with no subdirectory structure. This will be the absolute path prefix for the files listed in the manfiest
QCReads="/storage/home/stc5395/scratch/BordensteinVegVirome/qcreads" # this is optional
Manifest="/storage/home/stc5395/scratch/BordensteinVegVirome/VegViromeAssemblyManifest.tsv" # a tsv file with the following columns: SampleID, Lib1_F, Lib1_R, Lib2_F, Lib2_R where F and R are the forward and reverse reads, and lib1 and lib2 are if the sample was sequenced multiple times

#Example Manifest
#SampleID	Lib1_F	Lib1_R	Lib2_F	Lib2_R
#McReynolds1_C1_Control	C1_control_CKDN250003824-1A_22THGKLT4_L5_1.fq.gz	C1_control_CKDN250003824-1A_22THGKLT4_L5_2.fq.gz		
#McReynolds1_C2_Control	C2_control_CKDN250003825-1A_22TGMMLT4_L4_1.fq.gz	C2_control_CKDN250003825-1A_22TGMMLT4_L4_2.fq.gz	C2_control_CKDN250003825-1A_22THGKLT4_L5_1.fq.gz	C2_control_CKDN250003825-1A_22THGKLT4_L5_2.fq.gz
#McReynolds1_C3_Control	C3_control_CKDN250003826-1A_22THGHLT4_L3_1.fq.gz	C3_control_CKDN250003826-1A_22THGHLT4_L3_2.fq.gz	C3_control_CKDN250003826-1A_22THGKLT4_L5_1.fq.gz	C3_control_CKDN250003826-1A_22THGKLT4_L5_2.fq.gz


######################################################################################################################################################################################
# Do not modify below this line

#Set up directories and logs
mkdir -p $OutDir
mkdir -p $OutDir/logs

WorkDir=$(mktemp -d -p /scratch/${USER}) # putting temp working files on scratch to avoid storage overruns. The folder name is randomly generated
JobID=$SLURM_JOB_ID
ArrayID=$(expr $SLURM_ARRAY_TASK_ID + 1)
Nthreads=$SLURM_CPUS_PER_TASK 
MemReq=$((SLURM_MEM_PER_NODE / 1024)) # Returns memory values in MB, need to convert to GB, always returns an integer

#Test settings for interactive node
#JobID=1
#ArrayID=2
#Nthreads=4
#MemReq=32

SampleID=$(sed "${ArrayID}q;d" ${Manifest} | cut -f 1)
R1_L1=$(sed "${ArrayID}q;d" ${Manifest} | cut -f 2)
R2_L1=$(sed "${ArrayID}q;d" ${Manifest} | cut -f 3)
R1_L2=$(sed "${ArrayID}q;d" ${Manifest} | cut -f 4)
R2_L2=$(sed "${ArrayID}q;d" ${Manifest} | cut -f 5)
LibNumber=$(if [ "$R2_L2" = "" ]; then echo 1; else echo 2; fi) # to know if there are multiple libraries

echo "$(date) Running on $SampleID with $LibNumber lanes and files: $R1_L1 $R2_L1 $R1_L2 $R2_L2 outputing to $OutDir/${SampleID} and temporary storage in $WorkDir" >> $OutDir/logs/MGS_MasterLog.txt # dump starts to a master log
exec > ${OutDir}/logs/${SampleID}.log 2> ${OutDir}/logs/${SampleID}.err # move each sample to it's own log
echo "$(date) Running on $SampleID"
echo "Number of threads: $Nthreads" 
echo "Memory Usage: $MemReq"
echo "Number of Libraries: $LibNumber" 
echo "Reads: $R1_L1 $R2_L1 $R1_L2 $R2_L2" 
echo "Output directory:  $OutDir"
echo "Temporary directory: $WorkDir"
echo "Job ID: $JobID"
echo "Array ID: $ArrayID"


#############################################################################################################################################################################################
# Move reads to temp directory
mkdir -p $WorkDir/reads
cp $Reads/${R1_L1} $WorkDir/reads/${SampleID}_L1_R1.fastq.gz
cp $Reads/${R2_L1} $WorkDir/reads/${SampleID}_L1_R2.fastq.gz
if [ "$LibNumber" == "2" ]; then
	cp $Reads/${R1_L2} $WorkDir/reads/${SampleID}_L2_R1.fastq.gz
	cp $Reads/${R2_L2} $WorkDir/reads/${SampleID}_L2_R2.fastq.gz
fi

#############################################################################################################################################################################################
# QC via FastP v0.22.0
if $RunFastP; then
	echo $(date) Starting fastp
	for i in $(seq 1 $LibNumber); do
		mkdir -p $WorkDir/fastpout
		mkdir -p ${OutDir}/fastp
		mamba run -n fastp_env fastp \
		--in1 $WorkDir/reads/${SampleID}_L${i}_R1.fastq.gz \
		--in2 $WorkDir/reads/${SampleID}_L${i}_R2.fastq.gz \
		--out1 $WorkDir/fastpout/${SampleID}_L${i}_1.fastq.gz \
		--out2 $WorkDir/fastpout/${SampleID}_L${i}_2.fastq.gz \
		--detect_adapter_for_pe \
		--trim_poly_g \
		--cut_front \
		--cut_tail \
		--cut_window_size 4 \
		--cut_mean_quality 20 \
		--length_required 60 \
		--dedup \
		--json $OutDir/fastp/${SampleID}_L${i}.json \
		--html $OutDir/fastp/${SampleID}_L${i}.html \
		--thread $Nthreads
		
		mv $WorkDir/fastpout/${SampleID}_L${i}_1.fastq.gz $WorkDir/reads/${SampleID}_L${i}_R1.fastq.gz
		mv $WorkDir/fastpout/${SampleID}_L${i}_2.fastq.gz $WorkDir/reads/${SampleID}_L${i}_R2.fastq.gz
	done
fi


#############################################################################################################################################################################################
# Remove Host Reads using Bowtie2 v 2.5.2
# Note using human genome GCF_000001405.39_GRCh38.p13

if [ "$HostSpecies" == "mouse" ]; then 
	echo "Mouse genome not available at this time!"
	exit 1
elif [ "$HostSpecies" == "human" ]; then
	KNDB="/storage/home/stc5395/work/database/human_genome/HumanGenome"
fi

if $RemoveHost; then
	echo "$(date) Starting hostgenome removal"
	
	mkdir -p "$WorkDir/hosttrim"
	module load bowtie2
	
	for i in $(seq 1 $LibNumber); do
  	mkdir -p "$WorkDir/hosttrim"
  
    bowtie2 \
      --very-sensitive \
      --dovetail \
      -x "$KNDB" \
      -1 "$WorkDir/reads/${SampleID}_L${i}_R1.fastq.gz" \
      -2 "$WorkDir/reads/${SampleID}_L${i}_R2.fastq.gz" \
      --un-conc-gz "$WorkDir/hosttrim/${SampleID}_L${i}_paired_%.fastq.gz" \
      -p "$Nthreads" \
      -S "$WorkDir/hosttrim/${SampleID}_L${i}.sam"
    
    rm "$WorkDir/hosttrim/${SampleID}_L${i}.sam"
    
    mv "$WorkDir/hosttrim/${SampleID}_L${i}_paired_1.fastq.gz" \
       "$WorkDir/reads/${SampleID}_L${i}_R1.fastq.gz"
    
    mv "$WorkDir/hosttrim/${SampleID}_L${i}_paired_2.fastq.gz" \
       "$WorkDir/reads/${SampleID}_L${i}_R2.fastq.gz"
	 done
	 module unload bowtie2
fi

#############################################################################################################################################################################################
#Move Reads if needed
if $KeepQCReads; then
	mkdir -p $QCReads
	for i in $(seq 1 $LibNumber); do
		cp "$WorkDir/reads/${SampleID}_L${i}_R1.fastq.gz" "${QCReads}/${SampleID}_L${i}_R1.fastq.gz"
		cp "$WorkDir/reads/${SampleID}_L${i}_R2.fastq.gz" "${QCReads}/${SampleID}_L${i}_R2.fastq.gz"
	done
fi

#############################################################################################################################################################################################
# Merge reads if needed
if [ "$LibNumber" == "2" ]; then
	echo "$(date) Merging Reads"
	cat "$WorkDir/reads/${SampleID}_L1_R1.fastq.gz" \
	    "$WorkDir/reads/${SampleID}_L2_R1.fastq.gz" \
	    > "$WorkDir/reads/${SampleID}_R1.fastq.gz"

	cat "$WorkDir/reads/${SampleID}_L1_R2.fastq.gz" \
	    "$WorkDir/reads/${SampleID}_L2_R2.fastq.gz" \
	    > "$WorkDir/reads/${SampleID}_R2.fastq.gz"

else
	echo "$(date) Renaming Reads"
	mv "$WorkDir/reads/${SampleID}_L1_R1.fastq.gz" \
	   "$WorkDir/reads/${SampleID}_R1.fastq.gz"

	mv "$WorkDir/reads/${SampleID}_L1_R2.fastq.gz" \
	   "$WorkDir/reads/${SampleID}_R2.fastq.gz"
fi

echo "$(date) Merging Reads"
cat \
    "$WorkDir/reads/${SampleID}_R1.fastq.gz" \
    "$WorkDir/reads/${SampleID}_R2.fastq.gz" \
    > "$WorkDir/reads/${SampleID}_humann.fastq.gz"

#############################################################################################################################################################################################
# Microbecensus version 1.1.0
if $RunMicrobeCensus; then
	echo "$(date) Starting Average Genome Size Estimation"
	mkdir -p ${OutDir}/microbecensus
	mamba run -n microbecensus_env run_microbe_census.py -v -t $Nthreads $WorkDir/reads/${SampleID}_R1.fastq.gz,$WorkDir/reads/${SampleID}_R2.fastq.gz ${OutDir}/microbecensus/${SampleID}.microbecensus
fi

#############################################################################################################################################################################################
# Metaphlan v 4.2.2
# database is stored within the conda environment /storage/home/stc5395/work/.conda/envs/mpa/lib/python3.13/site-packages/metaphlan/metaphlan_databases

if $RunMetaphlan; then
	echo $(date) Starting metaphlan
	mkdir -p ${OutDir}/metaphlan
	
	mamba run -n mpa metaphlan \
	--output_file "$OutDir/metaphlan/${SampleID}.metaphlan" \
	--mapout "$WorkDir/${SampleID}.bowtie2.bz2" \
	--input_type fastq \
	--sample_id "$SampleID" \
	--nproc "$Nthreads" \
	--tmp_dir "$WorkDir/metaphlan_scratch" \
	-t rel_ab_w_read_stats \
	"$WorkDir/reads/${SampleID}_humann.fastq.gz"
fi

#############################################################################################################################################################################################
# Humann
if $RunHumann; then
	echo $(date) Starting Humann
	mkdir -p $WorkDir/humann
	humann \
	--input $WorkDir/reads/${SampleID}_humann.fastq.gz \
	--output $WorkDir/humann/$SampleID \
	--threads $Nthreads \
	--input-format fastq \
	--search-mode uniref90 \
	--verbose \
	--output-basename $SampleID
	
	mkdir -p $OutDir/humann
	cp $WorkDir/humann/${SampleID}/*tsv $OutDir/humann/
	cp $WorkDir/humann/${SampleID}/${SampleID}_humann_temp/${SampleID}_metaphlan_bugs_list.tsv $OutDir/humann/
	
fi

#############################################################################################################################################################################################
# Assemble with metaspades v 4.2.0
# needs own conda environment for python 3.8
if $RunMetaSpades; then
	echo $(date) Starting assembly
	mkdir -p "$WorkDir/assemblytmp"
	mkdir -p "$WorkDir/assemblies"
	mkdir -p "$OutDir/assemblies"
	if [ "$LibNumber" == "2" ]; then
		mamba run -n SPAdes_env /storage/home/stc5395/bin/SPAdes-4.2.0-Linux/bin/metaspades.py \
		--threads "$Nthreads" \
		--memory "$MemReq" \
		--tmp-dir "$WorkDir/assemblytmp" \
		--pe-1 1 "$WorkDir/reads/${SampleID}_L1_R1.fastq.gz" \
		--pe-2 1 "$WorkDir/reads/${SampleID}_L1_R2.fastq.gz" \
		--pe-1 2 "$WorkDir/reads/${SampleID}_L2_R1.fastq.gz" \
		--pe-2 2 "$WorkDir/reads/${SampleID}_L2_R2.fastq.gz" \
		-o "$WorkDir/assemblies/$SampleID"
	else
		mamba run -n SPAdes_env /storage/home/stc5395/bin/SPAdes-4.2.0-Linux/bin/metaspades.py \
		--threads "$Nthreads" \
		--memory "$MemReq" \
		--tmp-dir "$WorkDir/assemblytmp" \
		-1 "$WorkDir/reads/${SampleID}_R1.fastq.gz" \
		-2 "$WorkDir/reads/${SampleID}_R2.fastq.gz" \
		-o "$WorkDir/assemblies/$SampleID"
	fi
	cp "$WorkDir/assemblies/${SampleID}/scaffolds.fasta" "$OutDir/assemblies/${SampleID}_scaffolds.fasta"
	# Run QUAST
	mkdir -p "$OutDir/assemblies/assembly_stats"
	mkdir -p "$OutDir/assemblies/assembly_stats/${SampleID}"
	mamba run -n quast_env metaquast.py \
	--threads $Nthreads \
	--output-dir "$OutDir/assemblies/assembly_stats/${SampleID}" \
	--max-ref-number 0 \
	"$OutDir/assemblies/${SampleID}_scaffolds.fasta"
fi



#############################################################################################################################################################################################
echo $(date)" Finished $SampleID "
rm -r $WorkDir
```

## Binning, Dereplication, Taxonomic Classification, and Diversity Estimates

### Submission Script

```{bash, eval = F}
#!/usr/bin/env bash
# =============================================================================
# 00_submit_pipeline.sh — Master submission script
#
# Run this interactively from the login node:
#   bash 00_submit_pipeline.sh
#
# Pipeline:
#
#   01 (array) ──► 02 ──► 03 ──┬──► 04  GTDBtk
#                               ├──► 05 ──► 06  DRAM annotate → distill
#                               └──► 07 ──► 08 (array) ──► 09  MetaPop
#
# Steps 04, 05, and 07 all launch in parallel immediately after dRep (03).
# =============================================================================

SCRIPT_DIR="/storage/home/stc5395/scratch/BordensteinVegVirome/BinningTaxonomyMetapopScripts"

# Helper: submit and capture job ID
submit() {
  local label="$1"; shift
  local jid
  jid=$(sbatch --parsable "$@")
  echo "  Submitted ${label}: job ${jid}"
  echo "$jid"
}

echo "========================================"
echo "Submitting MAG pipeline"
echo "Working dir: $SCRIPT_DIR"
echo "========================================"

# ------------------------------------------------------------------
# Step 01: Per-sample array (bowtie2 → MetaBAT2 → CheckM2)
# ------------------------------------------------------------------
ARRAY_JID=$(submit "01_mag_pipeline (array)" \
  "${SCRIPT_DIR}/01_mag_pipeline.sh")

# ------------------------------------------------------------------
# Step 02: Stage bins + build genome_info.csv (after array)
# ------------------------------------------------------------------
PREP_JID=$(submit "02_prep_drep" \
  --dependency=afterok:${ARRAY_JID} \
  "${SCRIPT_DIR}/02_prep_drep.sh")

# ------------------------------------------------------------------
# Step 03: dRep dereplication (after staging)
# ------------------------------------------------------------------
DREP_JID=$(submit "03_drep_dereplicate" \
  --dependency=afterok:${PREP_JID} \
  "${SCRIPT_DIR}/03_drep_dereplicate.sh")

# ------------------------------------------------------------------
# Steps 04, 05, 07 run in PARALLEL after dRep
# ------------------------------------------------------------------

# Step 04: GTDBtk on representatives
GTDB_JID=$(submit "04_gtdbtk_classify" \
  --dependency=afterok:${DREP_JID} \
  "${SCRIPT_DIR}/04_gtdbtk_classify.sh")

# Step 05: DRAM annotate
# Note: if you want DRAM to embed GTDBtk taxonomy, add a manual
# dependency here: --dependency=afterok:${GTDB_JID}
# That will delay annotation until GTDBtk finishes.
DRAM_ANN_JID=$(submit "05_dram_annotate" \
  --dependency=afterok:${DREP_JID} \
  "${SCRIPT_DIR}/05_dram_annotate.sh")

# Step 07: Build combined MAG reference + bowtie2 index for MetaPop
METAPOP_PREP_JID=$(submit "07_metapop_prep" \
  --dependency=afterok:${DREP_JID} \
  "${SCRIPT_DIR}/07_metapop_prep.sh")

# ------------------------------------------------------------------
# Downstream steps (serial within each branch)
# ------------------------------------------------------------------

# Step 06: DRAM distill (after annotate)
DRAM_DISTILL_JID=$(submit "06_dram_distill" \
  --dependency=afterok:${DRAM_ANN_JID} \
  "${SCRIPT_DIR}/06_dram_distill.sh")

# Step 08: MetaPop remap array (after reference is built)
REMAP_JID=$(submit "08_metapop_remap (array)" \
  --dependency=afterok:${METAPOP_PREP_JID} \
  "${SCRIPT_DIR}/08_metapop_remap_array.sh")

# Step 09: MetaPop run (after all remaps complete)
METAPOP_JID=$(submit "09_metapop_run" \
  --dependency=afterok:${REMAP_JID} \
  "${SCRIPT_DIR}/09_metapop_run.sh")

echo ""
echo "========================================"
echo "All jobs submitted"
echo "========================================"
printf "  %-30s %s\n" "01_mag_pipeline (array)"   "$ARRAY_JID"
printf "  %-30s %s\n" "02_prep_drep"               "$PREP_JID"
printf "  %-30s %s\n" "03_drep_dereplicate"        "$DREP_JID"
printf "  %-30s %s\n" "04_gtdbtk_classify"         "$GTDB_JID"
printf "  %-30s %s\n" "05_dram_annotate"           "$DRAM_ANN_JID"
printf "  %-30s %s\n" "06_dram_distill"            "$DRAM_DISTILL_JID"
printf "  %-30s %s\n" "07_metapop_prep"            "$METAPOP_PREP_JID"
printf "  %-30s %s\n" "08_metapop_remap (array)"   "$REMAP_JID"
printf "  %-30s %s\n" "09_metapop_run"             "$METAPOP_JID"
echo ""
echo "Monitor with:  squeue -u \$USER"
echo "Cancel all:    scancel $ARRAY_JID $PREP_JID $DREP_JID $GTDB_JID $DRAM_ANN_JID $DRAM_DISTILL_JID $METAPOP_PREP_JID $REMAP_JID $METAPOP_JID"
```

### MAG Pipeline

```{bash, eval = F}
#!/usr/bin/env bash
#SBATCH --job-name=MAG_Pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=80GB
#SBATCH --array=1-523
#SBATCH --output=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/pipeline_%A_%a.out
#SBATCH --error=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/pipeline_%A_%a.err
#SBATCH --partition=sla-prio
#SBATCH --account=one_sc_default
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=stc5395@psu.edu

# =============================================================================
# Step 01: Per-sample bowtie2 → MetaBAT2 → CheckM2
# GTDBtk is intentionally excluded here; it runs on dereplicated
# representatives only (see 04_gtdbtk_classify.sh).
# =============================================================================

# ---------------------------
# USER CONFIGURABLES
# ---------------------------
BASE="/storage/home/stc5395/scratch/BordensteinVegVirome"
SAMPLE_LIST="${BASE}/samples.txt"
ASSEMBLY_DIR="${BASE}/assemblies"
BOWTIE2_OUT="${BASE}/bowtie2"
METABAT_OUT="${BASE}/metabat2"
CHECKM2_OUT="${BASE}/checkm2"
READ_DIR="${BASE}/qcreads"
# point to diamond db
CHECKM2_DB="/storage/home/stc5395/scratch/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"

METABAT_ENV="metabat2_env"
CHECKM2_ENV="checkm2"

# ---------------------------
# SOFTWARE
# ---------------------------
module load bowtie2
module load samtools

THREADS="${SLURM_CPUS_PER_TASK:-1}"

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d '\r')"
[[ -n "$SAMPLE" ]] || { echo "ERROR: empty sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

echo "==============================="
echo "SAMPLE: $SAMPLE"
echo "TASK:   ${SLURM_ARRAY_TASK_ID}"
echo "DATE:   $(date)"
echo "==============================="

# ============================================================
# 1) BOWTIE2 — map reads to own assembly (for MetaBAT2 depth)
# ============================================================
READ1="${READ_DIR}/${SAMPLE}_L1_R1.fastq.gz"
READ2="${READ_DIR}/${SAMPLE}_L1_R2.fastq.gz"
ASSEMBLY="${ASSEMBLY_DIR}/${SAMPLE}_scaffolds.fasta"

[[ -s "$READ1" ]]    || { echo "ERROR: missing READ1: $READ1"; exit 1; }
[[ -s "$READ2" ]]    || { echo "ERROR: missing READ2: $READ2"; exit 1; }
[[ -s "$ASSEMBLY" ]] || { echo "ERROR: missing ASSEMBLY: $ASSEMBLY"; exit 1; }

# ------------------------------------------------------------------
# Rename assembly contig headers: prepend sample ID to every sequence
# name so contigs are globally unique across all samples.
#
# SPAdes produces headers like:  >NODE_1_length_456_cov_12.34_ID_7
# After rename:                  >S001_NODE_1_length_456_cov_12.34_ID_7
#
# The renamed file is written alongside the original; the original is not
# modified. Skips the awk step if the file already exists.
# ------------------------------------------------------------------
ASSEMBLY_RENAMED="${ASSEMBLY_DIR}/rehead/${SAMPLE}_scaffolds.fasta"

if [[ ! -s "$ASSEMBLY_RENAMED" ]]; then
  echo "[$(date)] Renaming assembly headers for $SAMPLE"
  awk -v s="$SAMPLE" '/^>/{print ">" s "_" substr($0, 2); next} {print}' \
    "$ASSEMBLY" > "$ASSEMBLY_RENAMED"
  echo "Renamed: $(grep -c '^>' "$ASSEMBLY_RENAMED") contigs → $ASSEMBLY_RENAMED"
else
  echo "Renamed assembly already exists — skipping rename step"
fi

mkdir -p "${BOWTIE2_OUT}/${SAMPLE}"
cd "${BOWTIE2_OUT}/${SAMPLE}" || exit 1

LOG="${BOWTIE2_OUT}/${SAMPLE}/${SAMPLE}.mapping.log"
{
  echo "[$(date)] START $SAMPLE"
  echo "READ1=$READ1 | READ2=$READ2"
  echo "ASSEMBLY_RENAMED=$ASSEMBLY_RENAMED | THREADS=$THREADS"
  echo "[$(date)] bowtie2-build"
} > "$LOG"

bowtie2-build --threads "$THREADS" "$ASSEMBLY_RENAMED" "${SAMPLE}_idx" >>"$LOG" 2>&1
if [[ $? -ne 0 ]]; then
  echo "ERROR: bowtie2-build failed for $SAMPLE" >&2
  exit 1
fi

echo "[$(date)] bowtie2 | samtools sort" >>"$LOG"
bowtie2 --threads "$THREADS" -x "${SAMPLE}_idx" -1 "$READ1" -2 "$READ2" 2>>"$LOG" \
  | samtools view -b - 2>>"$LOG" \
  | samtools sort -@ "$THREADS" -o "${SAMPLE}.sorted.bam" 2>>"$LOG"
if [[ $? -ne 0 ]]; then echo "ERROR: bowtie2/sort failed for $SAMPLE" >&2; exit 1; fi

samtools index -@ "$THREADS" "${SAMPLE}.sorted.bam" >>"$LOG" 2>&1
if [[ $? -ne 0 ]]; then echo "ERROR: samtools index failed for $SAMPLE" >&2; exit 1; fi

echo "[$(date)] DONE mapping $SAMPLE" >>"$LOG"

# ============================================================
# 2) METABAT2
# ============================================================
BAM="${BOWTIE2_OUT}/${SAMPLE}/${SAMPLE}.sorted.bam"
[[ -s "$BAM" ]]              || { echo "ERROR: missing BAM: $BAM"; exit 1; }
[[ -s "$ASSEMBLY_RENAMED" ]] || { echo "ERROR: missing renamed FASTA: $ASSEMBLY_RENAMED"; exit 1; }

mkdir -p "${METABAT_OUT}/${SAMPLE}"
cd "${METABAT_OUT}/${SAMPLE}" || exit 1

MBLOG="${METABAT_OUT}/${SAMPLE}/${SAMPLE}.metabat2.log"
{
  echo "[$(date)] START $SAMPLE"
  echo "BAM=$BAM | FASTA=$ASSEMBLY_RENAMED | THREADS=$THREADS"
  echo "[$(date)] metabat2 initiation"
} > "$MBLOG"

mamba run -n "$METABAT_ENV" runMetaBat.sh \
  --minContig 2500 \
  --maxP 95 \
  --minS 60 \
  --minCV 1 \
  --minClsSize 200000 \
  --numThreads "$THREADS" \
  "$ASSEMBLY_RENAMED" \
  "$BAM" >>"$MBLOG" 2>&1

if [[ $? -ne 0 ]]; then
  echo "ERROR: MetaBAT2 failed for $SAMPLE" >&2
  echo "[$(date)] FAIL MetaBAT2" >>"$MBLOG"
  exit 1
fi

echo "[$(date)] DONE MetaBAT2 $SAMPLE" >>"$MBLOG"

# ============================================================
# HELPER: locate MetaBAT2 bins subdirectory
# runMetaBat.sh creates a timestamped dir: scaffolds.fasta.metabat-bins-XXXXXXXX/
# ============================================================
BINS_SUBDIR=$(find "${METABAT_OUT}/${SAMPLE}" \
  -maxdepth 1 -mindepth 1 -type d \
  \( -name "*.metabat-bins*" -o -name "bins" \) 2>/dev/null \
  | sort | tail -1)

if [[ -z "$BINS_SUBDIR" ]]; then
  echo "WARNING: no MetaBAT2 bins subdir found for $SAMPLE — falling back to parent dir"
  BINS_SUBDIR="${METABAT_OUT}/${SAMPLE}"
fi

BIN_COUNT=$(find "$BINS_SUBDIR" -maxdepth 1 -name "*.fa" 2>/dev/null | wc -l)
echo "Bins found: $BIN_COUNT in $BINS_SUBDIR"

# ============================================================
# 3) CHECKM2
# ============================================================
if [[ "$BIN_COUNT" -eq 0 ]]; then
  echo "WARNING: no bins for $SAMPLE — skipping CheckM2"
else
  mkdir -p "${CHECKM2_OUT}/${SAMPLE}"

  C2LOG="${CHECKM2_OUT}/${SAMPLE}/${SAMPLE}.checkm2.log"
  {
    echo "[$(date)] START $SAMPLE"
    echo "BINS_DIR=$BINS_SUBDIR ($BIN_COUNT bins) | THREADS=$THREADS"
    echo "[$(date)] checkm2 initiation"
  } > "$C2LOG"

  mamba run -n "$CHECKM2_ENV" checkm2 predict \
    --input "$BINS_SUBDIR" \
    --output-directory "${CHECKM2_OUT}/${SAMPLE}" \
    --force \
    --extension fa \
    --database_path "$CHECKM2_DB" \
    --threads "$THREADS" >>"$C2LOG" 2>&1

  if [[ $? -ne 0 ]]; then
    echo "ERROR: CheckM2 failed for $SAMPLE" >&2
    echo "[$(date)] FAIL CheckM2" >>"$C2LOG"
    exit 1
  fi

  echo "[$(date)] DONE CheckM2 $SAMPLE" >>"$C2LOG"
fi

echo "Pipeline complete for $SAMPLE at $(date)"
```

### Stage Bins

```{bash, eval = F}
#!/usr/bin/env bash
#SBATCH --job-name=Prep_dRep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --output=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/prep_drep_%j.out
#SBATCH --error=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/prep_drep_%j.err
#SBATCH --partition=sla-prio
#SBATCH --account=one_sc_default
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stc5395@psu.edu

# =============================================================================
# Step 02: Stage all MetaBAT2 bins with sample-prefixed names and build
# genome_info.csv from CheckM2 quality reports for dRep.
#
# Why rename bins:
#   Every sample produces bin.1.fa, bin.2.fa, etc.  Copying them into a
#   single flat directory without renaming causes collisions.  We use the
#   convention <SAMPLE>__<original_bin_name>.fa (double underscore as separator
#   so you can parse the sample ID back out later if needed).
#
# genome_info.csv:
#   dRep uses pre-computed CheckM2 quality scores to (a) filter out bins
#   below quality thresholds before dereplication and (b) select the highest-
#   quality genome as the representative for each cluster.
# =============================================================================

# ---------------------------
# USER CONFIGURABLES
# ---------------------------
BASE="/storage/home/stc5395/scratch/BordensteinVegVirome"
SAMPLE_LIST="${BASE}/samples.txt"
METABAT_OUT="${BASE}/metabat2"
CHECKM2_OUT="${BASE}/checkm2"
DREP_STAGING="${BASE}/drep/all_bins"
GENOME_INFO="${BASE}/drep/genome_info.csv"

# ---------------------------
mkdir -p "$DREP_STAGING" "$(dirname "$GENOME_INFO")" \
  /storage/home/stc5395/scratch/LCR/mag_pipeline/log

# dRep requires this exact CSV header
echo "genome,completeness,contamination" > "$GENOME_INFO"

TOTAL_BINS=0
SAMPLES_WITH_BINS=0
SAMPLES_NO_BINS=0
SAMPLES_MISSING_QC=0

while IFS= read -r SAMPLE || [[ -n "$SAMPLE" ]]; do
  SAMPLE=$(echo "$SAMPLE" | tr -d '\r\n')
  [[ -z "$SAMPLE" ]] && continue

  # Locate MetaBAT2 bins subdirectory (mirrors logic from 01_mag_pipeline.sh)
  BINS_SUBDIR=$(find "${METABAT_OUT}/${SAMPLE}" \
    -maxdepth 1 -mindepth 1 -type d \
    \( -name "*.metabat-bins*" -o -name "bins" \) 2>/dev/null \
    | sort | tail -1)

  if [[ -z "$BINS_SUBDIR" ]]; then
    echo "WARNING [$SAMPLE]: no MetaBAT2 bins subdirectory found — skipping"
    (( SAMPLES_NO_BINS++ ))
    continue
  fi

  QC_REPORT="${CHECKM2_OUT}/${SAMPLE}/quality_report.tsv"
  if [[ ! -f "$QC_REPORT" ]]; then
    echo "WARNING [$SAMPLE]: no CheckM2 quality_report.tsv at $QC_REPORT"
    echo "  Bins from this sample will be staged but absent from genome_info.csv."
    echo "  dRep will re-run CheckM internally for these (slower, uses CheckM1)."
    (( SAMPLES_MISSING_QC++ ))
  fi

  SAMPLE_BIN_COUNT=0

  while IFS= read -r BIN_FA; do
    BASENAME=$(basename "$BIN_FA")       # e.g. bin.1.fa
    BIN_ID="${BASENAME%.fa}"             # e.g. bin.1
    NEWNAME="${SAMPLE}__${BASENAME}"     # e.g. S001__bin.1.fa

    cp "$BIN_FA" "${DREP_STAGING}/${NEWNAME}"

    if [[ -f "$QC_REPORT" ]]; then
      # CheckM2 quality_report.tsv columns: Name(1) Completeness(2) Contamination(3) ...
      COMP=$(awk -F'\t' -v id="$BIN_ID" 'NR>1 && $1==id {print $2; exit}' "$QC_REPORT")
      CONT=$(awk -F'\t' -v id="$BIN_ID" 'NR>1 && $1==id {print $3; exit}' "$QC_REPORT")

      if [[ -n "$COMP" && -n "$CONT" ]]; then
        echo "${NEWNAME},${COMP},${CONT}" >> "$GENOME_INFO"
      else
        echo "WARNING [$SAMPLE]: no CheckM2 row for bin '$BIN_ID'"
      fi
    fi

    (( SAMPLE_BIN_COUNT++ ))
    (( TOTAL_BINS++ ))
  done < <(find "$BINS_SUBDIR" -maxdepth 1 -name "*.fa" 2>/dev/null | sort)

  if [[ "$SAMPLE_BIN_COUNT" -gt 0 ]]; then
    echo "[$SAMPLE]: staged $SAMPLE_BIN_COUNT bins"
    (( SAMPLES_WITH_BINS++ ))
  else
    echo "WARNING [$SAMPLE]: bins subdir exists but contains no .fa files"
    (( SAMPLES_NO_BINS++ ))
  fi

done < "$SAMPLE_LIST"

GENOME_INFO_COUNT=$(tail -n +2 "$GENOME_INFO" | wc -l)

echo ""
echo "============================================================"
echo "dRep Staging Summary"
echo "============================================================"
echo "Samples with bins staged   : $SAMPLES_WITH_BINS"
echo "Samples with no bins       : $SAMPLES_NO_BINS"
echo "Samples missing CheckM2    : $SAMPLES_MISSING_QC"
echo "Total bins staged          : $TOTAL_BINS"
echo "Bins in genome_info.csv    : $GENOME_INFO_COUNT"
echo "Staging dir                : $DREP_STAGING"
echo "Genome info CSV            : $GENOME_INFO"
echo "============================================================"

if [[ "$TOTAL_BINS" -eq 0 ]]; then
  echo "ERROR: no bins staged — check paths and that 01_mag_pipeline.sh completed"
  exit 1
fi
```

### dRep
```{bash, eval = F}
#!/usr/bin/env bash
#SBATCH --job-name=dRep_MAGs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120GB
#SBATCH --output=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/drep_%j.out
#SBATCH --error=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/drep_%j.err
#SBATCH --partition=sla-prio
#SBATCH --account=one_sc_default
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stc5395@psu.edu

# =============================================================================
# Step 03: dRep dereplication
#
# Algorithm:
#   1. Primary clustering at 90% ANI (MASH-based, fast) groups bins into
#      coarse clusters — any two bins in the same primary cluster are related
#      enough to warrant a detailed comparison.
#   2. Secondary clustering at 95% ANI (skani, accurate) within each
#      primary cluster defines species-level OTUs (sOTUs / MAG clusters).
#      95% is the standard species boundary in GTDB and the field broadly.
#   3. The highest-quality representative (scored by: completeness - 5×
#      contamination + 0.5×log(N50) + ...) is selected per secondary cluster.
#
# Quality filtering (-comp / -con) is applied BEFORE dereplication.
# Only bins passing thresholds enter the secondary clustering and are eligible
# to become representatives.  Bins failing QC are excluded entirely.
#
# Thresholds used here (adjust to your study's needs):
#   -comp 50   completeness ≥ 50%   (MQ MAG standard; raise to 90 for HQ only)
#   -con  10   contamination ≤ 10%
#   -pa   0.9  primary ANI = 90%
#   -sa   0.95 secondary ANI = 95% (species-level)
#   -nc   0.30 minimum genome coverage for a valid ANI comparison
# =============================================================================

# ---------------------------
# USER CONFIGURABLES
# ---------------------------
BASE="/storage/home/stc5395/scratch/BordensteinVegVirome"
DREP_STAGING="${BASE}/drep/all_bins"
DREP_OUT="${BASE}/drep/output"
GENOME_INFO="${BASE}/drep/genome_info.csv"

DREP_ENV="drep_env"   # conda env containing dRep

# ---------------------------
THREADS="${SLURM_CPUS_PER_TASK:-1}"

[[ -d "$DREP_STAGING" ]] || { echo "ERROR: staging dir not found: $DREP_STAGING"; exit 1; }
[[ -f "$GENOME_INFO" ]]  || { echo "ERROR: genome_info.csv not found: $GENOME_INFO"; exit 1; }

TOTAL_INPUT=$(find "$DREP_STAGING" -name "*.fa" | wc -l)
echo "[$(date)] Starting dRep on $TOTAL_INPUT bins"

mkdir -p "$DREP_OUT"

mamba run -n "$DREP_ENV" dRep dereplicate \
  -g "${DREP_STAGING}"/*.fa \
  --genomeInfo "$GENOME_INFO" \
  --S_algorithm skani \
  -pa 0.90 \
  -sa 0.95 \
  -nc 0.30 \
  -comp 50 \
  -con 10 \
  --multiround_primary_clustering \
  -p "$THREADS" \
  "$DREP_OUT" 2>&1

if [[ $? -ne 0 ]]; then
  echo "ERROR: dRep failed — check logs in ${DREP_OUT}/log/" >&2
  exit 1
fi

# ------------------------------------------------------------------
# Report
# ------------------------------------------------------------------
REP_COUNT=$(find "${DREP_OUT}/dereplicated_genomes" -name "*.fa" 2>/dev/null | wc -l)
echo ""
echo "============================================================"
echo "dRep Summary"
echo "============================================================"
echo "Input bins        : $TOTAL_INPUT"
echo "Representatives   : $REP_COUNT"
echo "Reduction         : $(echo "scale=1; (1 - $REP_COUNT / $TOTAL_INPUT) * 100" | bc)%"
echo "Representatives   : ${DREP_OUT}/dereplicated_genomes/"
echo "============================================================"

if [[ "$REP_COUNT" -eq 0 ]]; then
  echo "ERROR: no representative genomes produced — check dRep logs"
  exit 1
fi

echo "[$(date)] dRep complete"
```

### GTDB-tk
```{bash, eval = F}
#!/usr/bin/env bash
#SBATCH --job-name=GTDBtk_Reps
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=320GB
#SBATCH --output=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/gtdbtk_%j.out
#SBATCH --error=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/gtdbtk_%j.err
#SBATCH --partition=sla-prio
#SBATCH --account=one_sc_default
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stc5395@psu.edu

# =============================================================================
# Step 04: GTDBtk classify_wf on dRep dereplicated representatives
# =============================================================================

# ---------------------------
# USER CONFIGURABLES
# ---------------------------
BASE="/storage/home/stc5395/scratch/BordensteinVegVirome"
DREP_OUT="${BASE}/drep/output"
GTDBTK_OUT="${BASE}/gtdbtk"

GTDBTK_ENV="gtdbtk-2.7.2"

# Set this only if GTDBTK_DATA_PATH is not already exported in your conda
# activate script or .bashrc.  Leave blank otherwise.
GTDBTK_DATA_PATH=""

# ---------------------------
THREADS="${SLURM_CPUS_PER_TASK:-1}"
REPS_DIR="${DREP_OUT}/dereplicated_genomes"

[[ -d "$REPS_DIR" ]] || { echo "ERROR: dRep representatives dir not found: $REPS_DIR"; exit 1; }

REP_COUNT=$(find "$REPS_DIR" -name "*.fa" | wc -l)
if [[ "$REP_COUNT" -eq 0 ]]; then
  echo "ERROR: no .fa files found in $REPS_DIR"
  exit 1
fi
echo "[$(date)] Running GTDBtk on $REP_COUNT representative MAGs"

mkdir -p "$GTDBTK_OUT"

GTLOG="${GTDBTK_OUT}/gtdbtk.log"
{
  echo "[$(date)] START"
  echo "REPS_DIR=$REPS_DIR ($REP_COUNT genomes)"
  echo "THREADS=$THREADS"
} > "$GTLOG"

[[ -n "$GTDBTK_DATA_PATH" ]] && export GTDBTK_DATA_PATH

mamba run -n "$GTDBTK_ENV" gtdbtk classify_wf \
  --genome_dir "$REPS_DIR" \
  --out_dir "$GTDBTK_OUT" \
  --cpus "$THREADS" \
  --pplacer_cpus 4 \
  --extension fa \
  >>"$GTLOG" 2>&1

if [[ $? -ne 0 ]]; then
  echo "ERROR: GTDBtk failed — see $GTLOG" >&2
  exit 1
fi

echo "[$(date)] DONE GTDBtk" >>"$GTLOG"

# The primary output files are:
#   ${GTDBTK_OUT}/classify/gtdbtk.bac120.summary.tsv   (bacteria)
#   ${GTDBTK_OUT}/classify/gtdbtk.ar53.summary.tsv     (archaea)
echo ""
echo "GTDBtk complete. Taxonomy summaries:"
ls -lh "${GTDBTK_OUT}/classify/"*.summary.tsv 2>/dev/null || echo "(no summary files found — check log)"
```

### DRAM Annotate and Distill

```{bash, eval = F}
#!/usr/bin/env bash
#SBATCH --job-name=DRAM_Annotate
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=500GB
#SBATCH --output=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/dram_annotate_%j.out
#SBATCH --error=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/dram_annotate_%j.err
#SBATCH --partition=sla-prio
#SBATCH --account=one_sc_default
#SBATCH --time=120:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stc5395@psu.edu

# =============================================================================
# Step 05: DRAM annotate — functional annotation of dereplicated MAGs
#
# DRAM (Distilled and Refined Annotation of Metabolism) annotates MAGs
# against multiple databases simultaneously:
#   - KEGG (KO assignments, pathway membership)
#   - UniRef90 (broad functional homology)
#   - Pfam (domain annotation)
#   - dbCAN (CAZymes: carbohydrate-active enzymes)
#   - MEROPS (peptidase families)
#   - VOG (viral ortholog groups, relevant for prophages)
#
# The output annotations.tsv feeds directly into Step 06 (DRAM distill),
# which produces the cross-MAG metabolic summary.
#
# PREREQUISITES:
#   DRAM databases must be set up before running this script.
#   Check your configuration with: DRAM-setup.py print_config
#   If databases are not set up: DRAM-setup.py prepare_databases --output_dir /path/to/dram_dbs
#
# Memory note:
#   DRAM loads multiple large databases into memory simultaneously.
#   UniRef90 alone requires ~100 GB.  500 GB is a safe ceiling for a full
#   database suite.
#
# Time note:
#   Annotation runtime scales roughly linearly with the number of MAGs and
#   their total gene count.
# =============================================================================

# =============================================================================
# Step 06: DRAM distill — metabolic summary across all annotated MAGs
#
# The distill step takes the raw annotations.tsv from Step 05 and produces:
#
#   metabolism_summary.xlsx
#     The primary output.  An Excel workbook with one sheet per metabolic
#     module (carbon fixation, nitrogen cycling, sulfur cycling, etc.).
#     Each row is a gene family; each column is a MAG.  Presence/absence
#     and gene copy counts are summarized, making cross-MAG metabolic
#     comparison straightforward.
#
#   product.html
#     An interactive heatmap of metabolic modules across all MAGs, colored
#     by module completeness.  The most visually useful output for a first
#     pass at metabolic diversity.
#
#   genome_stats.tsv
#     Per-MAG summary statistics: genome size, GC%, coding density, rRNA
#     counts, tRNA counts, and completeness of key metabolic modules.
# =============================================================================

# ---------------------------
# USER CONFIGURABLES
# ---------------------------
BASE="/storage/home/stc5395/scratch/BordensteinVegVirome"
DREP_OUT="${BASE}/drep/output"
DRAM_ANNOTATE_OUT="${BASE}/dram/annotation"
DRAM_DISTILL_OUT="${BASE}/dram/distillation"

DRAM_ENV="DRAM_env"   # conda env containing DRAM v1.x

# Optional: path to GTDBtk taxonomy summary to embed taxonomy in DRAM output.
# Set to "" to skip.  GTDBtk must have completed (step 04) before using this.
GTDBTK_TAXONOMY="${BASE}/gtdbtk/classify/gtdbtk.bac120.summary.tsv"

# ---------------------------
THREADS="${SLURM_CPUS_PER_TASK:-1}"
REPS_DIR="${DREP_OUT}/dereplicated_genomes"

[[ -d "$REPS_DIR" ]] || { echo "ERROR: dRep representatives dir not found: $REPS_DIR"; exit 1; }

REP_COUNT=$(find "$REPS_DIR" -name "*.fa" | wc -l)
if [[ "$REP_COUNT" -eq 0 ]]; then
  echo "ERROR: no .fa files in $REPS_DIR"; exit 1
fi
echo "[$(date)] DRAM annotate on $REP_COUNT representative MAGs"

mkdir -p "$DRAM_ANNOTATE_OUT"

# Build optional taxonomy flag
TAXONOMY_FLAG=""
if [[ -n "$GTDBTK_TAXONOMY" && -f "$GTDBTK_TAXONOMY" ]]; then
  TAXONOMY_FLAG="--gtdb_taxonomy ${GTDBTK_TAXONOMY}"
  echo "Using GTDBtk taxonomy: $GTDBTK_TAXONOMY"
else
  echo "No GTDBtk taxonomy file provided or found — skipping taxonomy embedding"
fi

mamba run -n "$DRAM_ENV" DRAM.py annotate \
  -i "${REPS_DIR}/*.fa" \
  -o "$DRAM_ANNOTATE_OUT" \
  --threads "$THREADS" \
  --min_contig_size 2500 \
  $TAXONOMY_FLAG \
  2>&1

if [[ $? -ne 0 ]]; then
  echo "ERROR: DRAM annotate failed" >&2
  exit 1
fi

echo ""
echo "============================================================"
echo "DRAM Annotation complete"
echo "============================================================"
echo "Output dir      : $DRAM_ANNOTATE_OUT"
echo "Annotations     : ${DRAM_ANNOTATE_OUT}/annotations.tsv"
echo "rRNA calls      : ${DRAM_ANNOTATE_OUT}/rrnas.tsv"
echo "tRNA calls      : ${DRAM_ANNOTATE_OUT}/trnas.tsv"
echo ""
echo "Next: submit 06_dram_distill.sh"
echo "============================================================"
echo "[$(date)] DONE DRAM annotate"

# ---------------------------
[[ -f "${DRAM_ANNOTATE_OUT}/annotations.tsv" ]] || {
  echo "ERROR: annotations.tsv not found at ${DRAM_ANNOTATE_OUT}/annotations.tsv"
  echo "       Ensure 05_dram_annotate.sh completed successfully."
  exit 1
}

mkdir -p "$DRAM_DISTILL_OUT"

echo "[$(date)] Starting DRAM distill"

# Include rRNA and tRNA files if they were produced during annotation
RRNA_FLAG=""
TRNA_FLAG=""
[[ -f "${DRAM_ANNOTATE_OUT}/rrnas.tsv" ]] && RRNA_FLAG="--rrna_file ${DRAM_ANNOTATE_OUT}/rrnas.tsv"
[[ -f "${DRAM_ANNOTATE_OUT}/trnas.tsv" ]] && TRNA_FLAG="--trna_file ${DRAM_ANNOTATE_OUT}/trnas.tsv"

mamba run -n "$DRAM_ENV" DRAM.py distill \
  -i "${DRAM_ANNOTATE_OUT}/annotations.tsv" \
  -o "$DRAM_DISTILL_OUT" \
  $RRNA_FLAG \
  $TRNA_FLAG \
  2>&1

if [[ $? -ne 0 ]]; then
  echo "ERROR: DRAM distill failed" >&2
  exit 1
fi


echo ""
echo "============================================================"
echo "DRAM Distillation complete"
echo "============================================================"
echo "Output dir              : $DRAM_DISTILL_OUT"
echo "Metabolic summary       : ${DRAM_DISTILL_OUT}/metabolism_summary.xlsx"
echo "Interactive heatmap     : ${DRAM_DISTILL_OUT}/product.html"
echo "Genome stats            : ${DRAM_DISTILL_OUT}/genome_stats.tsv"
echo "============================================================"
echo "[$(date)] DONE DRAM distill"
```

### Metapop

```{bash, eval = F}
#!/usr/bin/env bash
#SBATCH --job-name=MetaPop_Prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --output=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/metapop_prep_%j.out
#SBATCH --error=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/metapop_prep_%j.err
#SBATCH --partition=sla-prio
#SBATCH --account=one_sc_default
#SBATCH --time=4:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=stc5395@psu.edu

# =============================================================================
# Step 07: Build the combined MAG reference for MetaPop
#
# MetaPop needs all samples mapped to the SAME reference genome set.
# That reference is the dereplicated representatives from dRep.
#
# This script:
#   1. Concatenates all representative .fa files into combined_mags.fasta,
#      renaming contig headers to <MAG_name>_<original_contig> to guarantee
#      uniqueness across MAGs (different MAGs can have contig_1, contig_2, etc.)
#   2. Builds the bowtie2 index used by the remap array (step 08)
# =============================================================================

# ---------------------------
# USER CONFIGURABLES
# ---------------------------
BASE="/storage/home/stc5395/scratch/BordensteinVegVirome"
DREP_OUT="${BASE}/drep/output"
METAPOP_BASE="${BASE}/metapop"
COMBINED_MAG_REF="${METAPOP_BASE}/combined_mags.fasta"
COMBINED_MAG_IDX="${METAPOP_BASE}/idx/combined_mags_idx"


# ---------------------------
THREADS="${SLURM_CPUS_PER_TASK:-1}"
REPS_DIR="${DREP_OUT}/dereplicated_genomes"

module load bowtie2

[[ -d "$REPS_DIR" ]] || { echo "ERROR: dRep representatives dir not found: $REPS_DIR"; exit 1; }
REP_COUNT=$(find "$REPS_DIR" -name "*.fa" | wc -l)
[[ "$REP_COUNT" -gt 0 ]] || { echo "ERROR: no .fa files found in $REPS_DIR"; exit 1; }
echo "[$(date)] Building combined reference from $REP_COUNT representative MAGs"

mkdir -p "${METAPOP_BASE}" "${METAPOP_BASE}/idx" \
  "${METAPOP_BASE}/bams" "${METAPOP_BASE}/norm" "${METAPOP_BASE}/output"

# ------------------------------------------------------------------
# 1) Concatenate representatives into a single FASTA
#    Rename headers: ><MAG_basename>_<original_contig_name>
#    This ensures uniqueness across all MAGs in the reference.
# ------------------------------------------------------------------
echo "[$(date)] Concatenating representatives → $COMBINED_MAG_REF"
> "$COMBINED_MAG_REF"   # truncate/create

while IFS= read -r FA; do
  MAG=$(basename "${FA%.fa}")
  awk -v mag="$MAG" '/^>/{print ">" mag "_" substr($0,2); next} {print}' "$FA"
done < <(find "$REPS_DIR" -name "*.fa" | sort) >> "$COMBINED_MAG_REF"

TOTAL_CONTIGS=$(grep -c '^>' "$COMBINED_MAG_REF")
echo "Combined reference written: $COMBINED_MAG_REF ($TOTAL_CONTIGS contigs)"

# ------------------------------------------------------------------
# 2) Build bowtie2 index
# ------------------------------------------------------------------
echo "[$(date)] Building bowtie2 index → $COMBINED_MAG_IDX"
bowtie2-build --threads "$THREADS" "$COMBINED_MAG_REF" "$COMBINED_MAG_IDX" 2>&1

if [[ $? -ne 0 ]]; then
  echo "ERROR: bowtie2-build failed" >&2
  exit 1
fi

echo ""
echo "============================================================"
echo "MetaPop reference preparation complete"
echo "============================================================"
echo "Combined FASTA  : $COMBINED_MAG_REF"
echo "bowtie2 index   : $COMBINED_MAG_IDX"
echo "MAGs in ref     : $REP_COUNT"
echo "Total contigs   : $TOTAL_CONTIGS"
echo "============================================================"
echo "[$(date)] DONE MetaPop prep"
```

```{bash, eval = F}
#!/usr/bin/env bash
#SBATCH --job-name=MetaPop_Remap
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80GB
#SBATCH --array=1-523
#SBATCH --output=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/metapop_remap_%A_%a.out
#SBATCH --error=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/metapop_remap_%A_%a.err
#SBATCH --partition=sla-prio
#SBATCH --account=one_sc_default
#SBATCH --time=12:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=stc5395@psu.edu

# =============================================================================
# Step 08: Remap all samples to the combined dRep MAG reference
#
# This is a SEPARATE re-mapping from Step 01.  Step 01 mapped each sample
# to its own SPAdes assembly (per-sample, for MetaBAT2 depth estimation).
# This step maps every sample to the SHARED dereplicated MAG reference so
# that MetaPop can compute cross-sample abundance and population diversity.
#
# BAM files are written to a flat directory (not per-sample subdirectories)
# because MetaPop's --input_samples expects all *.bam files in one directory.
# The filename without the .bam extension becomes the sample identifier in
# all MetaPop output tables and plots.
#
# NOTE ON SAM/BAM HANDLING:
# bowtie2 now writes an uncompressed .sam file directly (via -S), and the
# SAM -> BAM conversion, sorting, and indexing happen as separate, sequential
# samtools steps afterward. This makes each stage easier to debug and rerun
# independently, at the cost of extra disk I/O and temporary SAM disk usage.
# The intermediate SAM file is deleted after a successful sorted BAM is
# produced (set DELETE_SAM=false below to keep it).
# =============================================================================

# ---------------------------
# USER CONFIGURABLES
# ---------------------------
BASE="/storage/home/stc5395/scratch/BordensteinVegVirome"
SAMPLE_LIST="${BASE}/samples.txt"
READ_DIR="${BASE}/qcreads"
METAPOP_BASE="${BASE}/metapop"
COMBINED_MAG_REF="${METAPOP_BASE}/combined_mags.fasta"
COMBINED_MAG_IDX="${METAPOP_BASE}/idx/combined_mags_idx"
METAPOP_SAM_DIR="${METAPOP_BASE}/sams"
METAPOP_BAM_DIR="${METAPOP_BASE}/bams"
NORM_DIR="${METAPOP_BASE}/norm"

# Set to "false" to keep the intermediate .sam file after BAM creation
DELETE_SAM=true

# ---------------------------
module load bowtie2
module load samtools

THREADS="${SLURM_CPUS_PER_TASK:-1}"

SAMPLE="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$SAMPLE_LIST" | tr -d '\r')"
[[ -n "$SAMPLE" ]] || { echo "ERROR: empty sample for task ${SLURM_ARRAY_TASK_ID}"; exit 1; }

READ1="${READ_DIR}/${SAMPLE}_L1_R1.fastq.gz"
READ2="${READ_DIR}/${SAMPLE}_L1_R2.fastq.gz"

echo "MetaPop remap: $SAMPLE at $(date)"

[[ -s "$READ1" ]] || { echo "ERROR: missing READ1: $READ1"; exit 1; }
[[ -s "$READ2" ]] || { echo "ERROR: missing READ2: $READ2"; exit 1; }
[[ -s "$COMBINED_MAG_REF" ]] || {
  echo "ERROR: combined MAG reference not found: $COMBINED_MAG_REF"
  echo "       Ensure 07_metapop_prep.sh completed successfully."
  exit 1
}
[[ -f "${COMBINED_MAG_IDX}.1.bt2" ]] || {
  echo "ERROR: bowtie2 index not found (${COMBINED_MAG_IDX}.1.bt2)"
  echo "       Ensure 07_metapop_prep.sh completed successfully."
  exit 1
}

mkdir -p "$METAPOP_SAM_DIR" "$METAPOP_BAM_DIR" "$NORM_DIR"

SAM_FILE="${METAPOP_SAM_DIR}/${SAMPLE}.sam"
BAM_FILE_UNSORTED="${METAPOP_BAM_DIR}/${SAMPLE}.unsorted.bam"
BAM_FILE="${METAPOP_BAM_DIR}/${SAMPLE}.bam"

RLOG="${METAPOP_BAM_DIR}/${SAMPLE}.remap.log"
{
  echo "[$(date)] START $SAMPLE"
  echo "READ1=$READ1 | READ2=$READ2"
  echo "INDEX=$COMBINED_MAG_IDX | THREADS=$THREADS"
} > "$RLOG"

# ------------------------------------------------------------------
# Step 1: bowtie2 alignment -> SAM file
# ------------------------------------------------------------------
echo "[$(date)] bowtie2 alignment -> SAM" >>"$RLOG"

bowtie2 \
  --threads "$THREADS" \
  -x "$COMBINED_MAG_IDX" \
  -1 "$READ1" \
  -2 "$READ2" \
  -S "$SAM_FILE" \
  2>>"$RLOG"

if [[ $? -ne 0 || ! -s "$SAM_FILE" ]]; then
  echo "ERROR: bowtie2 alignment failed for $SAMPLE" >&2
  echo "[$(date)] FAIL bowtie2" >>"$RLOG"
  exit 1
fi

echo "[$(date)] bowtie2 complete: $SAM_FILE" >>"$RLOG"

# ------------------------------------------------------------------
# Step 2: SAM -> BAM conversion
# ------------------------------------------------------------------
echo "[$(date)] samtools view (SAM -> BAM)" >>"$RLOG"

samtools view -b -@ "$THREADS" -o "$BAM_FILE_UNSORTED" "$SAM_FILE" 2>>"$RLOG"

if [[ $? -ne 0 || ! -s "$BAM_FILE_UNSORTED" ]]; then
  echo "ERROR: SAM -> BAM conversion failed for $SAMPLE" >&2
  echo "[$(date)] FAIL samtools view" >>"$RLOG"
  exit 1
fi

# ------------------------------------------------------------------
# Step 3: sort BAM
# ------------------------------------------------------------------
echo "[$(date)] samtools sort" >>"$RLOG"

samtools sort -@ "$THREADS" -o "$BAM_FILE" "$BAM_FILE_UNSORTED" 2>>"$RLOG"

if [[ $? -ne 0 || ! -s "$BAM_FILE" ]]; then
  echo "ERROR: BAM sort failed for $SAMPLE" >&2
  echo "[$(date)] FAIL samtools sort" >>"$RLOG"
  exit 1
fi

rm -f "$BAM_FILE_UNSORTED"

# ------------------------------------------------------------------
# Step 4: index sorted BAM
# ------------------------------------------------------------------
echo "[$(date)] samtools index" >>"$RLOG"

samtools index -@ "$THREADS" "$BAM_FILE" 2>>"$RLOG" || {
  echo "ERROR: samtools index failed for $SAMPLE" >&2; exit 1
}

# ------------------------------------------------------------------
# Step 5: clean up intermediate SAM (optional)
# ------------------------------------------------------------------
if [[ "$DELETE_SAM" == "true" ]]; then
  rm -f "$SAM_FILE"
  echo "[$(date)] Removed intermediate SAM: $SAM_FILE" >>"$RLOG"
else
  echo "[$(date)] Retained intermediate SAM: $SAM_FILE" >>"$RLOG"
fi

# ------------------------------------------------------------------
# Norm file entry
# MetaPop --norm expects: <BAM_name_without_.bam>\t<total_reads>
# flagstat "in total" = all records (mapped + unmapped) = library size
# ------------------------------------------------------------------
TOTAL_READS=$(samtools flagstat "$BAM_FILE" \
  | awk '/in total/{print $1; exit}')
[[ -z "$TOTAL_READS" ]] && echo "WARNING: could not determine total reads for $SAMPLE"
echo -e "${SAMPLE}\t${TOTAL_READS}" > "${NORM_DIR}/${SAMPLE}.norm"
echo "[$(date)] Norm: $SAMPLE = $TOTAL_READS reads" >>"$RLOG"

echo "[$(date)] DONE remap $SAMPLE" >>"$RLOG"
echo "Done MetaPop remap $SAMPLE"
```

```{bash, eval = F}
#!/usr/bin/env bash
#SBATCH --job-name=MetaPop_Run
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=200GB
#SBATCH --output=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/metapop_%j.out
#SBATCH --error=/storage/home/stc5395/scratch/BordensteinVegVirome/logs/mag_pipeline/metapop_%j.err
#SBATCH --partition=sla-prio
#SBATCH --account=one_sc_default
#SBATCH --time=24:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=stc5395@psu.edu

# =============================================================================
# Step 09: Run MetaPop — population-level diversity analysis across MAGs
#
# MetaPop outputs:
#   Preprocessing:    per-MAG depth, breadth, and coverage across all samples
#   Macrodiversity:   normalized abundance, alpha/beta diversity across samples
#   Microdiversity:   SNP density, nucleotide diversity (pi), Watterson's
#                     theta, Tajima's D, pN/pS per MAG per sample
#   Visualization:    summary plots for all modules
#
# MetaPop filters MAGs with insufficient depth/breadth automatically;
# only MAGs detected at sufficient coverage in a sample contribute to
# microdiversity for that sample.
# =============================================================================

# ---------------------------
# USER CONFIGURABLES
# ---------------------------
BASE="/storage/home/stc5395/scratch/BordensteinVegVirome"
METAPOP_BASE="${BASE}/metapop"
COMBINED_MAG_REF="/storage/home/stc5395/scratch/BordensteinVegVirome/drep/output/dereplicated_genomes"
METAPOP_BAM_DIR="${METAPOP_BASE}/bams"
NORM_DIR="${METAPOP_BASE}/norm"
METAPOP_OUT="${METAPOP_BASE}/output"

METAPOP_ENV="metapop_env"

# MetaPop filtering thresholds
ID_MIN=95     # minimum % identity for read-to-reference alignment
MIN_COV=20    # minimum breadth of coverage (%) to include a MAG in microdiversity
MIN_DEP=10    # minimum truncated average depth for microdiversity

# ---------------------------
THREADS="${SLURM_CPUS_PER_TASK:-1}"

mkdir -p "$METAPOP_OUT"

# ------------------------------------------------------------------
# Combine per-sample norm files
# ------------------------------------------------------------------
NORM_FILE="${METAPOP_OUT}/norm.tsv"
echo "[$(date)] Combining per-sample norm files..."

NORM_COUNT=$(find "$NORM_DIR" -name "*.norm" 2>/dev/null | wc -l)
if [[ "$NORM_COUNT" -eq 0 ]]; then
  echo "ERROR: no .norm files in $NORM_DIR — did 08_metapop_remap_array.sh complete?"
  exit 1
fi
cat "${NORM_DIR}"/*.norm > "$NORM_FILE"
echo "Norm file: $NORM_FILE ($NORM_COUNT samples)"
awk -F'\t' '$2==0||$2=="" {print "WARNING: zero reads for sample:",$1}' "$NORM_FILE"

# ------------------------------------------------------------------
# Sanity checks
# ------------------------------------------------------------------
BAM_COUNT=$(find "$METAPOP_BAM_DIR" -name "*.bam" ! -name "*.bai" 2>/dev/null | wc -l)
echo "BAMs found: $BAM_COUNT"
[[ "$BAM_COUNT" -eq 0 ]] && { echo "ERROR: no BAMs in $METAPOP_BAM_DIR"; exit 1; }
[[ -s "$COMBINED_MAG_REF" ]] || { echo "ERROR: reference not found: $COMBINED_MAG_REF"; exit 1; }

TOTAL_CONTIGS=$(grep -c '^>' "$COMBINED_MAG_REF")
echo "[$(date)] Starting MetaPop: $BAM_COUNT samples, $TOTAL_CONTIGS reference contigs"

MPLOG="${METAPOP_OUT}/metapop.log"

mamba run -n "$METAPOP_ENV" metapop \
  --input_samples "$METAPOP_BAM_DIR" \
  --reference "$COMBINED_MAG_REF" \
  --norm "$NORM_FILE" \
  --threads "$THREADS" \
  --id_min "$ID_MIN" \
  --min_cov "$MIN_COV" \
  --min_dep "$MIN_DEP" \
  > "$MPLOG" 2>&1

if [[ $? -ne 0 ]]; then
  echo "ERROR: MetaPop failed — see $MPLOG" >&2
  exit 1
fi

echo ""
echo "============================================================"
echo "MetaPop complete"
echo "============================================================"
echo "Output : $METAPOP_OUT"
echo "Log    : $MPLOG"
echo "============================================================"
echo "[$(date)] DONE MetaPop"
```
