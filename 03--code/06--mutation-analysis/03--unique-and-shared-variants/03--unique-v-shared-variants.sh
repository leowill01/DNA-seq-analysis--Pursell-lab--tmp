#!/usr/bin/env zsh

# FILENAME:		03--unique-v-shared-variants.sh
# AUTHOR:		Leo Williams
# DATE:			Created: 2020 04_Apr 03
# USAGE:		[USAGE]
# DESCRIPTION:	See https://stackoverflow.com/questions/430078/shell-script-templates for template ideas
# VERSION:		[VERSION]
# NOTES:		[NOTES]
# ZSH VERSION:	[ZSH VERSION]
# DEV PLATFORM:	[DEV PLATFORM]

# ==============================================================================

# ARGUMENTS ########################################

INPUT=$1
OUTPUT_DIR=$2

# SETUP ########################################

# Make timestamped dir for results output
RESULTS="${OUTPUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir $RESULTS

# MAIN CODE ########################################

# bcftools isec --nfiles 2 *.vcf.gz | wc -l

for every **/*.vcf in input dir,
	- bgzip and tabix index;
	- 
