#!/usr/bin/env zsh

# FILENAME:		-
# AUTHOR:		Leo Williams
# DATE:			2020-04-22
# USAGE:		-
# DESCRIPTION:	See https://stackoverflow.com/questions/430078/shell-script-templates for template ideas
# VERSION:		-
# NOTES:		-
# ZSH VERSION:	-
# DEV PLATFORM:	macOS 10.15.4

# ==============================================================================

# ARGUMENTS ########################################

FWD_FQ=$1
REV_FQ=$2
OUT_PREFIX=$3
REF_FASTA=$4
OUTPUT_DIR=$5

# SETUP ########################################

# Make timestamped dir for results output
RESULTS_DIR="${OUTPUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir $RESULTS_DIR

# echo arguments for error log
{
	echo "GATK4 executable:" $GATK4_EXEC
	echo "FASTQ 1:" $FWD_FQ
	echo "FASTQ 2:" $REV_FQ
	echo "Output prefix:" $OUT_PREFIX
	echo "Reference genome FASTA:" $REF_FASTA
} &>>"${RESULTS_DIR}/error-log.txt"

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$RESULTS_DIR"

# set java version to 8
# jdk 1.8 # FIXME: need to source function

# MAIN CODE ########################################

# make uBAM from FASTQ files with FastqToSam
	# Only required options are shown
	# <SAMPLE_READ_GROUP_NAME> goes on uBAM read group header
	gatk FastqToSam \
	--FASTQ="$FWD_FQ" \
	--FASTQ2="$REV_FQ" \
	--OUTPUT="${RESULTS_DIR}/${OUT_PREFIX}.RG.u.bam" \
	--SAMPLE_NAME="$OUT_PREFIX" \
	--PLATFORM="illumina" \
	&>>"${RESULTS_DIR}/error-log.txt"
	# -TMP_DIR=$TMPDIR
	# --READ_GROUP_NAME="$OUT_PREFIX" \
	# --PLATFORM_UNIT="$OUT_PREFIX" \
	# --LIBRARY_NAME="$OUT_PREFIX" \

# rewrite BAM file with adapter-trimmed tags
	# This tool clears any existing adapter-trimming tags (XT:i:) in the optional tag region of a BAM file. The SAM/BAM file must be sorted by query name.
	# Outputs a metrics file histogram showing counts of bases_clipped per read.
	gatk MarkIlluminaAdapters \
	--INPUT="${RESULTS_DIR}/${OUT_PREFIX}.RG.u.bam" \
	--OUTPUT="${RESULTS_DIR}/${OUT_PREFIX}.RG.mkadap.u.bam" \
	--METRICS="${RESULTS_DIR}/${OUT_PREFIX}.mkadap.u.bam.metrics.txt" \
	&>>"${RESULTS_DIR}/error-log.txt"
	# -TMP_DIR=$TMPDIR

	# Remove intermediate input uBAM
	rm "${RESULTS_DIR}/${OUT_PREFIX}.RG.u.bam"

# convert Illumina adapter-marked uBAM back to FASTQ with SamToFastq in order to assign adapters low quality so they don't get mistaken for a read
	gatk SamToFastq \
	--INPUT="${RESULTS_DIR}/${OUT_PREFIX}.RG.mkadap.u.bam" \
	--FASTQ="${RESULTS_DIR}/${OUT_PREFIX}.mkadap.fq" \
	--CLIPPING_ATTRIBUTE=XT \
	--CLIPPING_ACTION=2 \
	--INTERLEAVE=true \
	--INCLUDE_NON_PF_READS=true \
	&>>"${RESULTS_DIR}/error-log.txt"
	# -TMP_DIR=$TMPDIR

	# NOTE: DO NOT Remove intermediate input uBAM because you have to merge it later with the aligned BAM

# align FASTQ sequence reads to a reference genome using BWA-MEM
	bwa mem -M -t 7 -p \
	"$REF_FASTA" \
	"${RESULTS_DIR}/${OUT_PREFIX}.mkadap.fq" > "${RESULTS_DIR}/${OUT_PREFIX}.bare.bam" \
	2>>"${RESULTS_DIR}/error-log.txt"

	# Remove intermediate adapter-marked FASTQ
	"${RESULTS_DIR}/${OUT_PREFIX}.mkadap.fq"

# Use GATK MergeBamAlignment to add metadata from uBAM back to the adapter-marked/removed/discounted and aligned BAM and to sort the final BAM by coordinate, which is required for use by MarkDuplicates
	gatk MergeBamAlignment \
	--REFERENCE_SEQUENCE="$REF_FASTA" \
	--UNMAPPED_BAM="${RESULTS_DIR}/${OUT_PREFIX}.RG.mkadap.u.bam" \
	--ALIGNED_BAM="${RESULTS_DIR}/${OUT_PREFIX}.bare.bam" \
	--OUTPUT="${RESULTS_DIR}/${OUT_PREFIX}.merged.sorted.bam" \
	--SORT_ORDER=coordinate \
	--CREATE_INDEX=true \
	--ADD_MATE_CIGAR=true \
	--CLIP_ADAPTERS=false \
	--CLIP_OVERLAPPING_READS=true \
	--INCLUDE_SECONDARY_ALIGNMENTS=true \
	--MAX_INSERTIONS_OR_DELETIONS=-1 \
	--PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
	--ATTRIBUTES_TO_RETAIN=XS \
	&>>"${RESULTS_DIR}/error-log.txt"
	# -TMP_DIR=$TMPDIR

	# Remove intermediate pre-merge uBAM and BAM files
	rm "${RESULTS_DIR}/${OUT_PREFIX}.RG.mkadap.u.bam" "${RESULTS_DIR}/${OUT_PREFIX}.bare.bam" "${RESULTS_DIR}/${OUT_PREFIX}.bare.bai"

# Mark read duplicates in BAM files that are artifacts from sequencing with MarkDuplicates (Picard).
	# It takes a previously-made adapter-removed, sorted, merged, and indexed BAM file as input.
	gatk MarkDuplicates \
	--INPUT="${RESULTS_DIR}/${OUT_PREFIX}.merged.sorted.bam" \
	--OUTPUT=results/"${RESULTS_DIR}/${OUT_PREFIX}.clean.bam" \
	--METRICS_FILE="${RESULTS_DIR}/${OUT_PREFIX}.clean.bam.metrics.txt" \
	--CREATE_INDEX=true \
	&>>"${RESULTS_DIR}/error-log.txt"
	# -TMP_DIR=$TMPDIR

	# Remove intermediate merged sorted BAM and .bai
	rm "${RESULTS_DIR}/${OUT_PREFIX}.merged.sorted.bam" "${RESULTS_DIR}/${OUT_PREFIX}.merged.sorted.bai"

echo "done"