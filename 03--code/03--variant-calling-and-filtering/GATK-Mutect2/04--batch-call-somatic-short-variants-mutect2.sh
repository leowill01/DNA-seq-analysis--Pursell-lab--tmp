#!/usr/bin/env bash

# title: Somatic variant calling in tumor-normal pairs with Mutect2
# author: Leo Williams
# date: 2020-09-18-17:18:54--Friday, Sep 18, 2020 05:18 PM

# USAGE ########################################

#

# ABOUT ########################################

# This script batch queues up the script to run Mutect2 from a 2-column table of tumor-normal paired BAMs to run the analyses on.

# INPUT:
# - 2-column TSV with absolute filepath to the tumor BAMs as the 1st column and absolute filepath to the normal BAMs in the 2nd column


# OUTPUT:
# - none

# ARGUMENTS ########################################

in_tumor_normal_table=$1 # required
in_mutect2_sbatch_script=$2 # required
in_ref_fasta=$3 # required
out_dir=$4 # required
# in_panel_of_normals_vcf_gz=$5

printf "Args:
in_tumor_normal_table: "$in_tumor_normal_table"\n\
in_mutect2_sbatch_script: "$in_mutect2_sbatch_script"\n\
in_ref_fasta: "$in_ref_fasta"\n\
out_dir: "$out_dir"\n"
# in_panel_of_normals_vcf_gz: "$in_panel_of_normals_vcf_gz""

# in_tumor_bam=$2 # required
# in_tumor_bam_basename=$(basename "$in_tumor_bam")
# in_normal_bam=$3 # required
# in_normal_bam_basename=$(basename "$in_normal_bam")
# normal_sample_name=$4 # required
# out_prefix="T-${in_tumor_bam_basename%.bam}--N-${in_normal_bam_basename%.bam}"
# in_panel_of_normals_vcf_gz=$6 # optional
# in_germline_resource=$7 # optional

# SETUP ########################################

# Make timestamped dir for results output
dir_results="${out_dir}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir "$dir_results"

# Make a copy of the script run and put it into the results dir
path_script=$(echo "$0")
cp "$path_script" "$dir_results"

# MAIN CODE ########################################

## BATCH QUEUE SCRIPT
# remove 1st line of in_tumor_normal_table
headerless_tn_table="${dir_results}/headerless-tumor-normal-table.tsv"
sed 1d "$in_tumor_normal_table" > "$headerless_tn_table"


while read sample_row_line; do
	table_tumor_bam_path="$(echo "$sample_row_line" | cut -f 2)"
	table_tumor_bam_RG_sample_name="$(echo "$sample_row_line" | cut -f 3)"
	table_normal_bam_path="$(echo "$sample_row_line" | cut -f 4)"
	table_normal_bam_RG_sample_name="$(echo "$sample_row_line" | cut -f 5)"

	# make run-specific out_dir
	out_dir_run="${dir_results}/${table_tumor_bam_RG_sample_name}-${table_normal_bam_RG_sample_name}"
	mkdir "$out_dir_run"

	echo "Batch queueing Mutect2 without Panel of Normals..."

	sbatch "$in_mutect2_sbatch_script" \
	"$in_ref_fasta" \
	"$table_tumor_bam_path" \
	"$table_normal_bam_path" \
	"$table_normal_bam_RG_sample_name" \
	"$out_dir_run"
done < "$headerless_tn_table"

	# FIXME: disabled to see if problem was with PON arg
	# if PON is supplied, run with
	# if [[ -z "$in_panel_of_normals_vcf_gz" ]] ; then
	# 	echo "Batch queueing Mutect2 with Panel of Normals..."

	# 	sbatch "$in_mutect2_sbatch_script" \
	# 	"$in_ref_fasta" \
	# 	"$table_tumor_bam_path" \
	# 	"$table_normal_bam_path" \
	# 	"$table_normal_bam_RG_sample_name" \
	# 	"$out_dir_run" \
	# 	"$in_panel_of_normals_vcf_gz"
	# else
	# 	echo "Batch queueing Mutect2 without Panel of Normals..."

	# 	sbatch "$in_mutect2_sbatch_script" \
	# 	"$in_ref_fasta" \
	# 	"$table_tumor_bam_path" \
	# 	"$table_normal_bam_path" \
	# 	"$table_normal_bam_RG_sample_name" \
	# 	"$out_dir_run"
	# fi
# done < "$in_tumor_normal_table"
