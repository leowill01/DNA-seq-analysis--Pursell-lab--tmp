#!/usr/bin/env bash

# title: Batch queue Mutect2 tumor-only mode for PON creation
# author: Leo Williams
# date: 2020-09-18-17:18:54--Friday, Sep 18, 2020 05:18 PM

# USAGE ########################################

#

# ABOUT ########################################

# This script is part of the GATK pipeline for calling somatic variants with Mutect2. This script is part 1 of 2 steps in the larger first step of creating a Panel of Normals (PON) to use in Mutect2 to help mask artifactual variants present in the cohort.

# INPUT:
# - reference genome FASTA
# - list normal BAMs for PON
# - path to Mutect2 tumor-only script

# OUTPUT:
# - 

# ARGUMENTS ########################################

in_ref_fasta=$1
in_list_of_normal_bams=$2
in_mutect_tumor_only_script=$3
out_dir=$4
# out_pon_rootname=$5

# SETUP ########################################

# Make timestamped dir for results output
dir_results="${out_dir}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir "$dir_results"

# make subdir for GenomicsDB
dir_genomicsDB="${dir_results}/GenomicsDB"
mkdir "$dir_genomicsDB"

# Make a copy of the script run and put it into the results dir
path_script=$(echo "$0")
cp "$path_script" "$dir_results"

# MAIN CODE ########################################

## LOAD MODULES
module load gatk/4.1.8.1

## RUN ANALYSIS

# 1. Take each normal BAM in the list and queue it to run Mutect2 in tumor-only mode
while read normal_bam; do
	# sbatch "$in_mutect_tumor_only_script" \
	# "$in_ref_fasta" `# arg 1: in_ref_fasta` \
	# "$normal_bam" `# arg 2: in_sample_normal_bam` \
	# "$dir_results" `# arg 3: out_dir`
	sbatch "$in_mutect_tumor_only_script" \
	"$in_ref_fasta" \
	"$normal_bam" \
	"$dir_results"
done < "$in_list_of_normal_bams"

# END ########################################
echo "End"
