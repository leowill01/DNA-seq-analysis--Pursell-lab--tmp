#!/usr/bin/env zsh

# title: 05--filter-mutect-for-PASS.sh
# author: Leo Williams
# date: 2020-10-08

# USAGE ########################################

# ./script.sh $dirh_with_mutect_filtered_vcfs

# ABOUT ########################################

#

# ARGUMENTS ########################################

in_dirh_of_mutect_filtered_vcfs=$1
out_dir=$2

# SETUP ########################################

# Make timestamped dir for results output
dir_results="${out_dir}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir $dir_results

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$path_script" "$dir_results"

# MAIN CODE ########################################

for vcf in "$in_dirh_of_mutect_filtered_vcfs"/**/*.mutect2.filtered*.vcf; do
	vcf_basename=$(basename "$vcf")
	vcf_rootname="${vcf_basename%.vcf}"
	out_vcf_basename="${vcf_rootname}.filtPass.vcf"

	bcftools view -f '.,PASS' "$vcf" > "$dir_results"/"$out_vcf_basename"
done
