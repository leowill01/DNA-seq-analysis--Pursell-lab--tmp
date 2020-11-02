#!/usr/bin/env bash

# About ----------------------------------------

# This script scans a directory hierarchy for all VCF files and runs each individually through a MutationalPatterns.r script that analyzes a single VCF using the MutationalPatterns R package.
# This may require curating a folder of specific VCFs files to be analyzed in one run (but still individually). E.g. if you wanted to run multiple VCFs from the same sample but from different levels of filtering, or if you wanted to run 

# Usage ====

# bash batch-mut-sigs.sh [INPUT_VCF_DIR]

# Within the sample folder, make an `INPUT/` folder that contains the VCFs you want to individually analyze

# Setup ----------------------------------------

# Arguments ====

INPUT_VCF_DIR=$1
PATTERN=$2

# Execution ----------------------------------------

# for every VCF file in the INPUT_VCF_DIR hierarchy, 
for vcf in ${INPUT_VCF_DIR}/**/${PATTERN}; do
#     make the unique "vcf-name/" dir
    echo $(basename $vcf)
    mkdir $(basename $vcf)
#     change to the unique "vcf-name/" as the working dir
    cd $(basename $vcf)
#     run the mut-sigs.r script from within each vcf-specific working dir
    Rscript "../01--mutation-signatures-single.r" $vcf $(basename $vcf)
    # Change back to top dir
    cd ..
done