#!/bin/zsh
# Author: Leo Williams
# Date: 2019-11-26
# Platform: macOS Catalina v10.15.1

# Usage: 

# About ########################################

	# What: This script takes as input a two-column newline-separated list of amino acid changes in the 1st column and TCGA case ID barcodes in the 2nd column and outputs a tab-separated list of just the TCGA case IDs.

	# Why: 

	# How: 

# Setup ########################################

	# Change working dir to project root dir
	cd ../../..

	# Variables (script)

	# Variables (interactive)
	VAR_AND_TCGA_ID_LIST="02--data/01--raw-data/01--internal-Pursell/TCGA-ID-list-Hodel-CRISPR-paper/patient-sample-list.txt"
	SAMPLE_DIR="/Volumes/Pursell-Lab-HD/Pursell-lab/projects/project_02--DNA-seq-analysis-POLE-mutant-tumors/04--analysis/group--Hodel-2019-paper-TCGA-POLE-tumors/expt-00--obtain-TCGA-data"

	# Make output dir
	OUTPUT_DIR="${SAMPLE_DIR}/output--$(date +\%Y-\%m-\%d-%H%M%S)"
	mkdir $OUTPUT_DIR
	# All output files should have format `"${OUTPUT_DIR}/filename.txt"`. Use `basename` 

# Execution Code ########################################

	# read the variant-TCGA barcode list

	# extract just the case IDs in the 2nd column

	

# Session Information ########################################

