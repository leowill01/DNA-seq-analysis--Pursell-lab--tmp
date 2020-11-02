#!/bin/bash
#SBATCH --job-name=map-clean	                ### Job name
#SBATCH --output=Output_Log-map-clean.log	### File to store output
#SBATCH --error=Error_Log-map-clean.err	    ### File to store error messages
#SBATCH --qos=normal				    ### Quality of service queue
#SBATCH --time=24:00:00					### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1						### Nodes requested for job
#SBATCH --ntasks-per-node=1				### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1				### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G						### Memory requested for job. Max 64GB/node
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwilli24@tulane.edu

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# This script maps FASTQ reads to a reference genome and cleans up the aligned BAM according to the GATK4 best practices for DNAseq analysis.
# Tools used: GATK4/Picard Tools, BWA

# TODO: make results validation criteria and dir
# FIXME: does the bam .bai index appear in the "results/" dir?


################################################################################
# ARGUMENTS
################################################################################
GATK4_EXEC="/lustre/project/zpursell/leo/code-packages-external/gatk-4.1.2.0/gatk"
FWD_FQ=$1
REV_FQ=$2
OUT_PREFIX=$3
REF_FASTA=$4

echo "GATK4 executable:" $GATK4_EXEC
echo "FASTQ 1:" $FWD_FQ
echo "FASTQ 2:" $REV_FQ
echo "Output prefix:" $OUT_PREFIX
echo "Reference genome FASTA:" $REF_FASTA

################################################################################
# EXECUTION CODE
################################################################################
# Load Cypress modules
module load java-openjdk/1.8.0 # load Java version required for GATK4
module load bwa

# Make results dir
mkdir -p results

########################################
# Use GATK4 FastqToSam to make a uBAM from FASTQ files
# Only required options are shown
# <SAMPLE_READ_GROUP_NAME> goes on uBAM read group header
$GATK4_EXEC FastqToSam \
--FASTQ="$FWD_FQ" \
--FASTQ2="$REV_FQ" \
--OUTPUT="${OUT_PREFIX}.RG.u.bam" \
--SAMPLE_NAME="$OUT_PREFIX" \
--READ_GROUP_NAME="$OUT_PREFIX" \
--PLATFORM="illumina" \
-TMP_DIR=$TMPDIR
# --PLATFORM_UNIT="$OUT_PREFIX" \
# --LIBRARY_NAME="$OUT_PREFIX" \

########################################
# Read BAM file and rewrite it with new adapter-trimming tags.
# This tool clears any existing adapter-trimming tags (XT:i:) in the optional tag region of a BAM file. The SAM/BAM file must be sorted by query name.
# Outputs a metrics file histogram showing counts of bases_clipped per read.
$GATK4_EXEC MarkIlluminaAdapters \
--INPUT="${OUT_PREFIX}.RG.u.bam" \
--OUTPUT="${OUT_PREFIX}.RG.mkadap.u.bam" \
--METRICS="results/${OUT_PREFIX}.mkadap.u.bam.metrics.txt" \
-TMP_DIR=$TMPDIR

# Remove intermediate input uBAM
rm "${OUT_PREFIX}.RG.u.bam"

########################################
# Uses GATK4 SamToFastq to convert Illumina adapter-marked uBAM back to FASTQ in order to assign adapters low quality so they don't get mistaken for a read
$GATK4_EXEC SamToFastq \
--INPUT="$OUT_PREFIX".RG.mkadap.u.bam \
--FASTQ="$OUT_PREFIX".mkadap.fq \
--CLIPPING_ATTRIBUTE=XT \
--CLIPPING_ACTION=2 \
--INTERLEAVE=true \
--INCLUDE_NON_PF_READS=true \
-TMP_DIR=$TMPDIR

# NOTE: DO NOT Remove intermediate input uBAM because you have to merge it later with the aligned BAM

########################################
# Aligns FASTQ sequence reads to a reference genome using BWA-MEM
bwa mem -M -t 7 -p \
"$REF_FASTA" \
"$OUT_PREFIX".mkadap.fq > "$OUT_PREFIX".bare.bam

# Remove intermediate adapter-marked FASTQ
rm "$OUT_PREFIX".mkadap.fq

########################################
# Use GATK MergeBamAlignment to add metadata from uBAM back to the adapter-marked/removed/discounted and aligned BAM and to sort the final BAM by coordinate, which is required for use by MarkDuplicates
$GATK4_EXEC MergeBamAlignment \
--REFERENCE_SEQUENCE="$REF_FASTA" \
--UNMAPPED_BAM="$OUT_PREFIX".RG.mkadap.u.bam \
--ALIGNED_BAM="$OUT_PREFIX".bare.bam \
--OUTPUT="$OUT_PREFIX".merged.sorted.bam \
--SORT_ORDER=coordinate \
--CREATE_INDEX=true \
--ADD_MATE_CIGAR=true \
--CLIP_ADAPTERS=false \
--CLIP_OVERLAPPING_READS=true \
--INCLUDE_SECONDARY_ALIGNMENTS=true \
--MAX_INSERTIONS_OR_DELETIONS=-1 \
--PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
--ATTRIBUTES_TO_RETAIN=XS \
-TMP_DIR=$TMPDIR

# Remove intermediate pre-merge uBAM and BAM files
rm "$OUT_PREFIX".RG.mkadap.u.bam "$OUT_PREFIX".bare.bam "$OUT_PREFIX".bare.bai

########################################
# Mark read duplicates in BAM files that are artifacts from sequencing with MarkDuplicates (Picard).
# It takes a previously-made adapter-removed, sorted, merged, and indexed BAM file as input.
$GATK4_EXEC MarkDuplicates \
--INPUT="$OUT_PREFIX".merged.sorted.bam \
--OUTPUT=results/"$OUT_PREFIX".clean.bam \
--METRICS_FILE="results/${OUT_PREFIX}.clean.bam.metrics.txt" \
--CREATE_INDEX=true \
-TMP_DIR=$TMPDIR

# Remove intermediate merged sorted BAM and .bai
rm "$OUT_PREFIX".merged.sorted.bam
rm "$OUT_PREFIX".merged.sorted.bai

################################################################################
# END JOB
################################################################################
echo End Job
