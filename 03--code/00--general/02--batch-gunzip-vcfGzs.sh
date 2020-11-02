#!/bin/zsh

# About ########################################

    # This script takes a directory as input, recursively finds and unzips all *.vcf.gz files within the hierarchy, and outputs them to a new directory.

# Setup ########################################

    # args
    INPUT_DIR_OF_VCFS=$1
    OUTPUT_DIR=$2

    # # Make timestamped run results dir inside output dir
	RESULTS_DIR="${OUTPUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
	mkdir $RESULTS_DIR

# Analysis ########################################

    # Recursively find all *.vcf.gz files and gunzip them to the output dir
        # for every *vcf.gz file in the hierarchy, 
        for i in "$INPUT_DIR_OF_VCFS"/**/*.vcf.gz
        do
            # echo VCF filename
            echo "$i"
            basename "$i" '.gz'
            # gunzip to *.vcf while keeping original file and output to output dir
            gunzip -c "$i" > "${RESULTS_DIR}/$(basename ${i} '.gz')"
        done
