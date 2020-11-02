#!/usr/bin/env bash

# title: 01.01--batch-sbatch-VarscanSomaticMutations.sh
# author: Leo Williams
# date: 2020-07-07 6:17 PM

# USAGE ########################################

# bash 01.01--batch-sbatch-VarscanSomaticMutations.sh \
# "$IN_DIR_OF_SAMPLE_DIRS" \
# "$OUT_DIR_FOR_SAMPLES" \
# "$SBATCH_SCRIPT" \
# "$NORMAL_BAM"

# ABOUT ########################################

# run from anywhere, but IN_DIR_OF_SAMPLE_DIRS must have this structure:

# IN_DIR_OF_SAMPLE_DIRS/
# |-- sample1/
# |   |-- results--YYYY-MM-DD/
# |       |-- sample1.clean.bam
# |-- sample2/
# |   |-- results--YYYY-MM-DD/
# |       |-- sample2.clean.bam
# |-- script--batch-sbatch.sh
# |-- script--sbatch.sh

# This script runs all samples against a single normal defined by $NORMAL_BAM

# ARGUMENTS ########################################

IN_DIR_OF_SAMPLE_DIRS=$1
OUT_DIR_FOR_SAMPLES=$2
SBATCH_SCRIPT=$3
NORMAL_BAM=$4
SBATCH_SCRIPT_BASE=$(basename "$SBATCH_SCRIPT")
REF_FASTA="/lustre/project/zpursell/leo/projects/project_02--DNA-seq-analysis/02--data/02--reference-data/genome/mm10--UCSC/mm10_UCSC.fa"

printf "Arguments: \n \
IN_DIR_OF_SAMPLE_DIRS: ${IN_DIR_OF_SAMPLE_DIRS} \n \
OUT_DIR_FOR_SAMPLES: ${OUT_DIR_FOR_SAMPLES} \n \
SBATCH_SCRIPT: ${SBATCH_SCRIPT} \n \
SBATCH_SCRIPT_BASE: ${SBATCH_SCRIPT_BASE} \n \
REF_FASTA: ${REF_FASTA} \n \
NORMAL_BAM: ${NORMAL_BAM}"

# SETUP ########################################

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$SBATCH_SCRIPT" "$OUT_DIR_FOR_SAMPLES"

# MAIN CODE ########################################

# for each sample dir, schedule an sbatch script for that sample
for i_IN_SAMPLE_DIR in "$IN_DIR_OF_SAMPLE_DIRS"/*/ ; do
	# dirname without trailing '/'
	IN_SAMPLE_DIR_PATH=$(echo "${i_IN_SAMPLE_DIR%/}") # gives absolute path without trailing '/' if in dir is given as arg instead of sourcing from that dir

	# get sample name
	SAMPLE_NAME="$(basename "$i_IN_SAMPLE_DIR")"

	# make dir for sample in master out dir
	OUT_SAMPLE_DIR="$OUT_DIR_FOR_SAMPLES"/"$SAMPLE_NAME"
	mkdir "$OUT_SAMPLE_DIR"

	# copy sbatch script to sample out dir
	cp "$SBATCH_SCRIPT" "$OUT_SAMPLE_DIR"

	# cd to sample out dir
	cd "$OUT_SAMPLE_DIR"

	# printf "Loop args:\n \
	# IN_SAMPLE_DIR_PATH: ${IN_SAMPLE_DIR_PATH}\n \
	# SAMPLE_NAME: ${SAMPLE_NAME}\n \
	# OUT_SAMPLE_DIR: ${OUT_SAMPLE_DIR}"

	# que up the sbatch script with the relevant arguments
		# sbatch script args:
			# VARSCAN2_JARFILE="/lustre/project/zpursell/leo/code-packages-external/varscan-2.4.3/VarScan.v2.4.3.jar"
			# TUMOR_BAM=$1
			# NORMAL_BAM=$2 # defined as NORMAL_BAM
			# OUT_PREFIX=$3
			# OUTPUT_DIR=$4
			# REF_GENOME_FASTA=$5 # defined as REF_FASTA

		# Log the command used to run sbatch script
		printf "\n"
		printf "\nCommand entered to run sbatch script:\nsbatch ${SBATCH_SCRIPT_BASE} \\ \n ${IN_SAMPLE_DIR_PATH}/results/*.bam \\ \n ${NORMAL_BAM} \\ \n ${SAMPLE_NAME} \\ \n ${OUT_SAMPLE_DIR} \\ \n ${REF_FASTA}\n"

		# run sbatch script for this sample/iteration
		sbatch "$SBATCH_SCRIPT_BASE" \
		"$IN_SAMPLE_DIR_PATH"/results/*.bam \
		"$NORMAL_BAM" \
		"$SAMPLE_NAME" \
		"$OUT_SAMPLE_DIR" \
		"$REF_FASTA"
done
