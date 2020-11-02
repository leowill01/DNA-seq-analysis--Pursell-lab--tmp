#!/usr/bin/env bash

# run from anywhere, but IN_DIR_OF_SAMPLE_DIRS must have this structure:
#IN_DIR_OF_SAMPLE_DIRS/
#|-- sample1/
#|   |-- fwd_1.fq
#|   |-- rev_2.fq
#|-- sample2/
#|   |-- fwd_1.fq
#|   |-- rev_2.fq
#|-- batch-sbatch.sh
#|-- sbatch.sh

IN_DIR_OF_SAMPLE_DIRS=$1
OUT_DIR=$2
SBATCH_SCRIPT=$3
SBATCH_SCRIPT_BASE=$(basename "$SBATCH_SCRIPT")
REF_FASTA=$4

for iIN_SAMPLE_DIR in "$IN_DIR_OF_SAMPLE_DIRS"/*/ ; do
	# dirname without trailing '/'
	IN_SAMPLE_DIR_PATH=$(echo "${iIN_SAMPLE_DIR%/}") # gives absolute path without trailing '/' if in dir is given as arg instead of sourcing from that dir

	# get sample name
	SAMPLE_NAME="$(basename "$iIN_SAMPLE_DIR")"

	# make dir for sample in master out dir
	OUT_SAMPLE_DIR="$OUT_DIR"/"$SAMPLE_NAME"
	mkdir "$OUT_SAMPLE_DIR"

	# copy sbatch script to sample out dir
	cp "$SBATCH_SCRIPT" "$OUT_SAMPLE_DIR"

	# FIXME: does this run each script in its respective folder correctly?
	# cd to sample out dir
	cd "$OUT_SAMPLE_DIR"

	# que up the sbatch script with the relevant arguments
	# args:
	# FWD_FQ=$1
	# REV_FQ=$2
	# OUT_PREFIX=$3
	# REF_FASTA=$4
	sbatch "$SBATCH_SCRIPT_BASE" \
	"$IN_SAMPLE_DIR_PATH"/*-1fwd.fq \
	"$IN_SAMPLE_DIR_PATH"/*-2rev.fq \
	"$SAMPLE_NAME" \
	"$REF_FASTA"
done
