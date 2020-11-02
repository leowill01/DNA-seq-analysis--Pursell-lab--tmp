#!/bin/bash

# About ##############################

# This script runs an R script which analyzes a set of VCFs for mutation spectra patterns using the MutationalPatterns R package.

# Usage ##############################

# bash mut-sigs.sh <in_VCF_dir>



# Arguments ##############################

IN_VCF_TOP_DIR=$1

# Code ##############################

# Enable globstar to search for VCF files within input dir
shopt -s globstar # FIXME: not working!

# For a directory with multiple VCF files recursively within the hierarchy (e.g. the "../results/" dir output from expt_04-variant-annotation)

# For every annotated VCF file in that hierarchy
# TODO: implement globstar
for i in $(find $(echo $(pwd)) -name *.vcf); do
	printf "${i}\n"
done

# Run 01--mutation-signatures.R on it by passing the VCF file as an command line arg to R

# Save results in "sample01/vcfBasename/results"
