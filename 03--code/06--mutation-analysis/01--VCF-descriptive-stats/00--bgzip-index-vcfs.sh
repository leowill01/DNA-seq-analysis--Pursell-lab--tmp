#!/usr/bin/env zsh

# ABOUT ########################################

# In: a dir of .vcf files
# Out: .vcf.gz files with corresponding index files

# SETUP ########################################

IN_DIR_OF_VCFS=$1
OUT_DIR_GZ=$2


# ANALYSIS ########################################

for i in "$IN_DIR_OF_VCFS"/**/*.vcf
do
    VCF_BASENAME=$(basename "$i")
    GZ_FILE="${OUT_DIR_GZ}/${VCF_BASENAME}.gz"
    bgzip -c "$i" > "$GZ_FILE"
    tabix "$GZ_FILE"
done


