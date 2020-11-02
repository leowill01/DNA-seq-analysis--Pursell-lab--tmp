#!/usr/bin/env zsh

# title: 00--format-Varscan-VCF-DP4-to-typeInt.sh
# author: Leo Williams
# date: 2020-08-22

# USAGE ########################################

# script.sh <IN_DIR_OF_VCFS> <OUT_DIR>

# ABOUT ########################################

# This script takes as input VCFs generated from Varscan tumor-normal somatic mutation calling and reformats the FORMAT/DP4 field to be of type integer with 4 values instead of type string

# ARGUMENTS ########################################

IN_DIR_OF_VCFS=$1 # VCFs must be output by Varscan2
# IN_DIR_OF_VCFS="/Volumes/Pursell-Lab-HD/Pursell-lab/02--projects/project_02--DNA-seq-analysis/04--analysis/03--variant-calling--postprocessing--annotation/BRCA1-related/2020-06-10--BRCA1-KO-batch2-targetseq/VarScan2/01.01--combined-somatic-muts-only/00--annotated-vcfs-only-all-runs/02--vcfs-isec"
OUT_DIR=$2
# OUT_DIR="/Volumes/Pursell-Lab-HD/Pursell-lab/02--projects/project_02--DNA-seq-analysis/04--analysis/03--variant-calling--postprocessing--annotation/BRCA1-related/2020-06-10--BRCA1-KO-batch2-targetseq/VarScan2/01.01--combined-somatic-muts-only/00--annotated-vcfs-only-all-runs/03--format-Varscan-VCF-header-DP4-to-integer"

# SETUP ########################################

# Make timestamped dir for results output
RESULTS_DIR="${OUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir $RESULTS_DIR

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$RESULTS_DIR"

# MAIN CODE ########################################
# for each VCF in input dir, use sed to change the FORMAT ID for DP4 to be of type Integer with 4 fields (changed from type String with 1 field as default set by Varscan)
for vcf in $IN_DIR_OF_VCFS/**/*.vcf ; do
	echo "$vcf" # print filename
	vcfBASENAME=$(basename "$vcf") # get basename & remove path

	# assign new output filepath
	vcfNew="${RESULTS_DIR}/${vcfBASENAME%.vcf}.fmtDP4ID.vcf"

	# replace ID with sed
	sed 's/##FORMAT=<ID=DP4,Number=1,Type=String/##FORMAT=<ID=DP4,Number=4,Type=Integer/' "$vcf" > "$vcfNew"
done
