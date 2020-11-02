#!/bin/bash

# run from anywhere, but argument passed as input dir or files must have this structure:
# ./
# |---dir1/
# |   |---file1
# |---dir2/
#     |---file2

# Args of script that is being batch-queued:
# ARG1=
# ARG2=

IN_DIR_OF_SAMPLE_DIRS=$1
OUT_DIR_OF_SAMPLE_DIRS=$2
SBATCH_SCRIPT=$3
SBATCH_SCRIPT_BASE=$(basename "$SBATCH_SCRIPT")
REF_FASTA=$4

# args must satisfy the args of the sbatch script
# ARG1=""

# run sbatch script for every sample
for SAMPLE_DIR_NAME in "$IN_DIR_OF_SAMPLE_DIRS"/*/ ; do
	# dirname without trailing '/'
	IN_SAMPLE_DIR=$(echo "${SAMPLE_DIR_NAME%/}") # gives absolute path without trailing '/' if in dir is given as arg instead of sourcing from that dir

	# get sample name
	SAMPLE_NAME="$(basename "$IN_SAMPLE_DIR")"

	# make output dir for sample in master out dir
	OUT_SAMPLE_DIR="$OUT_DIR_OF_SAMPLE_DIRS"/"$SAMPLE_NAME"
	mkdir "$OUT_SAMPLE_DIR"

	# copy sbatch script to sample out dir for review purposes
	cp "$SBATCH_SCRIPT" "$OUT_SAMPLE_DIR"

	# cd to and run the sbatch script from the sample out dir to output results folder in correct location (rather than whatever dir you run this script from)
	cd "$OUT_SAMPLE_DIR"

	# from within each sample folder, run the sbatch script with the correct args as if you were running it individually
	sbatch "$SBATCH_SCRIPT_BASE" \
	"$IN_SAMPLE_DIR"/"$ARG1" \
	"$ARG2"
done
