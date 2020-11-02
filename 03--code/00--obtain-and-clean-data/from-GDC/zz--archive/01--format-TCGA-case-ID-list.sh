#!/bin/zsh
# Author: Leo Williams
# Date: 2019-11-06
# Platform: macOS

# About ########################################

	# This script takes as input a newline-separated list of The Cancer Genome Atlas (TCGA) barcodes which have the 'sample' number at the end and:

	# 1. Removes the sample number, giving just the barcode up to the participant ID
	# 2. Converts the newline-separated case ID list to a TSV, which is able to be input into TCGA as a case set to further filter

# Usage ########################################

	# Run this script interactively with the working dir as the project root dir

# Setup ########################################

	# Change working dir to project dir
	cd /Volumes/Pursell-Lab-HD/Pursell-lab/projects/project_05--REV1-tumors

	# Assign variables (script)
		# TCGA case ID list
		CASE_ID_LIST=$1
		# output dir
		OUTPUT_SAMPLE_DIR=$2
	# # Assign variables (interactive)
	# 	# TCGA case ID list
	# 	CASE_ID_LIST="02--data/TCGA-REV1-tumor-expression-IDs/00.01--REV1-bot10pct-underexpressed-subjID.txt"
	# 	OUTPUT_SAMPLE_DIR="04--analysis/01--data-procurement/01--format-TCGA-case-ID-list"

	# Make results dir
		RESULTS_DIR="${OUTPUT_SAMPLE_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
		mkdir $RESULTS_DIR
		# All output files should have format `"${RESULTS_DIR}/filename.txt"`. Use `basename` 

# Code ########################################

	# Read in the case ID list
	cat $CASE_ID_LIST |
	# Remove trailing TCGA sample ID from the rest of the case ID
	sed 's/\-[0-9]\{2\}$//' |
	# Replace newlines with tabs to produce a TSV required for loading into TCGA query
	paste -s - > "${RESULTS_DIR}/${$(basename $CASE_ID_LIST)%.*}.no-sample.tsv"
	# pastes with tab as delimiter by default

# Session info ########################################

	# Export script to results dir?
	cp 03--code/01--data-preprocessing/01--format-TCGA-case-ID-list.sh \
	$RESULTS_DIR
