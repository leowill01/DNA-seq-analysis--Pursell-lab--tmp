#!/bin/zsh

# Usage: zsh script.sh [in_dir] [file_ext] [out_dir]

# About ########################################

	# This script takes as input a dir of VCF files downloaded from the Genomic Data Commons (GDC) and:

	# 1. Outputs a list of their relative filepaths (to be used downstream)
	# 2. Extracts the GDC VCF UUID from the VCF filename to be used to match against a table to extract the GDC case/patient ID to be used as the sample name

#  Setup ########################################

	# variables
	IN_DIR=$1
	FILE_EXT=$2
	OUT_DIR=$3

	# Make results dir
	RESULTS_DIR="${OUT_DIR}/results--filebaserootnames--$(date +\%Y-\%m-\%d-%H%M%S)"
	mkdir $RESULTS_DIR


# Code ########################################

	# define function to get absolute filepath
	absfilepath () {
        echo "$(cd "$(dirname "$1")"; pwd -P)/$(basename "$1")"
    }

	# init empty file to store filepaths
	FILEPATHS_FILE="${RESULTS_DIR}/abs-filepaths.tsv"
	# add column header to file
	echo "absolute_filepath" > "$FILEPATHS_FILE"
	# loop through *."$FILE_EXT" files and add their absolute filepaths to the file
	for i in "$IN_DIR"/**/*."$FILE_EXT"; do
		absfilepath "$i" >> "$FILEPATHS_FILE"
	done

	# init empty file to store basenames
	BASENAME_FILE="${RESULTS_DIR}/file-basenames.tsv"
	# init empty file to store rootnames
	ROOTNAME_FILE="${RESULTS_DIR}/file-rootnames.tsv"
	# add column header to basename file
	echo "file_basename" > "$BASENAME_FILE"
	# add column header to rootname file
	echo "file_rootname" > "$ROOTNAME_FILE"
	# loop through *."$FILE_EXT" files and add their basenames and rootnames to the respective files
	for i in "$IN_DIR"/**/*."$FILE_EXT"; do
		BASENAME=$(basename "$i")
		echo "$BASENAME" >> "$BASENAME_FILE"
		ROOTNAME="${BASENAME%%.*}"
		echo "$ROOTNAME" >> "$ROOTNAME_FILE"
	done

	# combine both files into a TSV of both columns
	## init empty combined file
	COMBINED_FILE="${RESULTS_DIR}/all-filebaserootnames.tsv"
	# paste each line in succession from each file into the combined file
	paste "$FILEPATHS_FILE" "$BASENAME_FILE" "$ROOTNAME_FILE"> "$COMBINED_FILE"
