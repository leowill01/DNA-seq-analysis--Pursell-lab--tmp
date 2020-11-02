#!/bin/zsh
# Author: Leo Williams
# Date: 2019-12-02
# Platform: macOS

# Usage: 

# About ########################################

	# What: 

	# Why: 

	# How: 

# Setup ########################################

    # Variables ====
    INPUT_VCF_DIR=$1
    OUTPUT_DIR=$2
   
   	# Make output dir
	RESULTS_DIR="${OUTPUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
	mkdir $RESULTS_DIR
	# All output files should have format `"${RESULTS_DIR}/filename.txt"`. Use `basename` 

# Execution Code ########################################

    # extract only SNVs from the VCFs and save them as new VCFs
    SNP_DIR="${RESULTS_DIR}/snp"
    mkdir "$SNP_DIR"
    for i in "$INPUT_VCF_DIR"/**/*.vcf; do
        BASENAME=$(basename "$i" '.vcf')
        bcftools view \
        --types snps \
        "$i" > "${SNP_DIR}/${BASENAME}.snp.vcf"
    done

    # # confirm no indels are in the SNV files
    # SNV_CHECK="${RESULTS_DIR}/snp-indel-check.txt"
    # echo '' > "$SNV_CHECK"
    # for i in **/*.snp.vcf
    # do
    #     echo "$i" >> "$SNV_CHECK"
    #     bcftools stats "$i" | grep -iE -e 'number of indels:' >> "$SNV_CHECK"
    # done

    # extract only INDELs from the VCFs and save them as new VCFs
    INDEL_DIR="${RESULTS_DIR}/indel"
    mkdir "$INDEL_DIR"
    for i in "$INPUT_VCF_DIR"/**/*.vcf; do
        BASENAME=$(basename "$i" '.vcf')
        bcftools view \
        --types indels \
        "$i" > "${INDEL_DIR}/${BASENAME}.indel.vcf"
    done

    # # confirm no SNVs are in the Indel files
    # INDEL_CHECK="${RESULTS_DIR}/indel-snp-check.txt"
    # echo '' > "$INDEL_CHECK"
    # for i in **/*.snp.vcf
    # do
    #     echo "$i" >> "$INDEL_CHECK"
    #     bcftools stats "$i" | grep -iE -e 'number of snps:' >> "$INDEL_CHECK"
    # done


# Session Information ########################################
