#!/bin/bash

# run from anywhere, but argument passed as input dir or files must have this structure:
# ./
# |---dir1/
# |   |---file1
# |---dir2/
#     |---file2

# Args of script that is being batch-queued:
# REF_FASTA=$1
# IN_BAM=$2
# OUTPUT_DIR=$3
# OUT_PREFIX=$4
# VARSCAN2_JARFILE="/lustre/project/zpursell/leo/code-packages-external/varscan-2.4.3/VarScan.v2.4.3.jar"

IN_DIR_OF_SAMPLE_DIRS=$1
OUT_DIR_OF_SAMPLE_DIRS=$2
SBATCH_SCRIPT=$3
SBATCH_SCRIPT_BASE=$(basename "$SBATCH_SCRIPT")
REF_FASTA=$4

# args must satisfy the args of the sbatch script
# ARG1=""

# run sbatch script for every sample
for BAM in "$IN_DIR_OF_SAMPLE_DIRS"/*/results/*.bam ; do
	# print filename
	echo "FILENAME: ${BAM}"

	# make basename
	BASENAME=$(basename "$BAM")
	echo "BASENAME: ${BASENAME}"
	ROOTNAME_SAMPLENAME=$(echo "${BASENAME%%.*}")
	echo "ROOTNAME: ${ROOTNAME_SAMPLENAME}"
	echo ""

	# make output dir for sample in master out dir
	OUT_SAMPLE_DIR="$OUT_DIR_OF_SAMPLE_DIRS"/"$ROOTNAME_SAMPLENAME"
	mkdir "$OUT_SAMPLE_DIR"

	# copy sbatch script to sample out dir for review purposes
	cp "$SBATCH_SCRIPT" "$OUT_SAMPLE_DIR"

	# cd to and run the sbatch script from the sample out dir to output results folder in correct location (rather than whatever dir you run this script from)
	cd "$OUT_SAMPLE_DIR"

	# from within each sample folder, run the sbatch script with the correct args as if you were running it individually
	# script args:
		# REF_FASTA=$1
		# IN_BAM=$2
		# OUTPUT_DIR=$3
		# OUT_PREFIX=$4
		# VARSCAN2_JARFILE="/lustre/project/zpursell/leo/code-packages-external/varscan-2.4.3/VarScan.v2.4.3.jar"

		# queue script
		sbatch "$SBATCH_SCRIPT_BASE" \
		"$REF_FASTA" \
		"$BAM" \
		"$OUT_SAMPLE_DIR" \
		"$ROOTNAME_SAMPLENAME"
done
