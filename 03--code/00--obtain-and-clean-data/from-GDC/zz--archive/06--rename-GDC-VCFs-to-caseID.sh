#!/usr/bin/env zsh

# FILENAME:		[FILENAME]
# AUTHOR:		[AUTHOR]
# DATE:			[DATE]
# VERSION:		[VERSION]
# NOTES:		[NOTES]
# ZSH VERSION:	[ZSH VERSION]
# DEV PLATFORM:	[DEV PLATFORM]

# ==============================================================================

# USAGE ########################################

# zsh 06--rename-GDC-VCFs-to-caseID.sh $

# ABOUT ########################################

# This script takes input of a table with full filepath, basename, and rootnames of a set of GDC VCFs, as well as columns for matched GDC case IDs to the VCF file UUIDs (the file rootnames). It renames the VCFs with the case ID as the file rootname for each file. This makes it easier to plot case IDs as samples in mutation analysis.

# ARGUMENTS ########################################

IN_DIR_OF_GDC_VCFS=$1
OUT_DIR=$2
# IN_TAB_FILEBASEROOTNAMES_CASEIDUUID=

# SETUP ########################################

# Make timestamped dir for results output
RESULTS_DIR="${OUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir $RESULTS_DIR

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$RESULTS_DIR"

# MAIN CODE ########################################

# rename every .vcf file basename to case ID
for i in "${IN_DIR_OF_GDC_VCFS}"/**/*.vep.vcf ; do
    echo "Working on: $(basename "$i")"

    # grep case id from vcf file
    CASE_ID=$(egrep -m 1 -i "individual" "$i" | \
    egrep -io 'TCGA-[[:alnum:]]{2}-[[:alnum:]]{4}')

    # make new filename from case ID found in file
    NEW_FILENAME="${OUT_DIR}/${CASE_ID}.vep.vcf"

    # rename file
    mv "$i" "$NEW_FILENAME"

done

