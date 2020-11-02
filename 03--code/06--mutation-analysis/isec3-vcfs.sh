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

# ABOUT ########################################

# This script takes as input a folder of VCFs for multiple samples, but where there are 3 (multiple? - generalized) iterations for each sample denoted by the suffix on the VCF basenames. E.g.:
# - sample1-run1.vcf
# - sample1-run2.vcf
# - sample1-run3.vcf
# are 3 iterations of the same sample.

# On a per-sample basis, this script takes the intersection of the first two iterations, then takes the intersection of that result with the third iteration VCF and outputs a final VCF that is the intersection of all 3 VCFs but which only has the sample FORMAT information from the first sample.

# ARGUMENTS ########################################

IN_DIR_OF_SAMPLE_DIRS=$1 # with each sample dir containing iterated VCFs. iterated VCF rootnames must end in '1', '2', or '3'.
# OUT_DIR=$2
FILE_ITER_SUFFIX_DELIM="--" # with number behind it, eg 'sample--1.vcf', where '--' would be the suffix delimiter

# SETUP ########################################

# Make timestamped dir for results output
# RESULTS_DIR="${OUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
# mkdir $RESULTS_DIR

# Make a copy of the script run and put it into the results dir
# SCRIPT_PATH=$(echo "$0")
# cp "$SCRIPT_PATH" "$RESULTS_DIR"

# MAIN CODE ########################################

# bgzip and tabix all VCFs in hierarchy
for i in "$IN_DIR_OF_SAMPLE_DIRS"/**/*.vcf; do
	bgzip -c "$i" > "${i}.gz"
	tabix "${i}.gz"
done

# this could all be a one-off to organize the files by sample dirs
	# # make a list of sample names and make a directory for each
	# for i in *.vcf ; do
	# 	echo "${i%-*}"
	# done | uniq | xargs mkdir

	# # for every sample directory, move all files with sample name to within
	# for isample_dir in * ; do
	# 	if [[ -d "$isample_dir" ]] ; then
	# 		# echo "$isample_dir"
	# 		# assign sample name
	# 		isample_name="$isample_dir"

	# 		# for every VCF with that sample name, move to the dir
	# 		for ivcf in "$isample_name"*.vcf ; do
	# 			mv "$ivcf" "$isample_dir"
	# 		done
	# 	fi
	# done
#

# for each sample dir with 3 VCF files, take 2 intersections and output a final 3-isec VCF with variants only common to all 3 files
for iSAMPLE_DIR in "$IN_DIR_OF_SAMPLE_DIRS"/* ; do
	# echo "$iSAMPLE_DIR"

	# cd to sample dir
	cd "$iSAMPLE_DIR"
	printf \
	"\n \
	Current directory: " ; pwd
	printf "\n"

	# assign files as variables
	FILE1=$(echo "$iSAMPLE_DIR"/*"$FILE_ITER_SUFFIX_DELIM"1*.vcf.gz)
	FILE2=$(echo "$iSAMPLE_DIR"/*"$FILE_ITER_SUFFIX_DELIM"2*.vcf.gz)
	FILE3=$(echo "$iSAMPLE_DIR"/*"$FILE_ITER_SUFFIX_DELIM"3*.vcf.gz)

	# print sample dir and files:
	printf \
	"Sample Dir: ${iSAMPLE_DIR}\n\
	File 1: ${FILE1}\n\
	File 2: ${FILE2}\n\
	File 3: ${FILE3}\n\n"

	# take isec of first two files
	bcftools isec \
	-p "$iSAMPLE_DIR"/"01--isec1-2" \
	"$FILE1" \
	"$FILE2"

	# bgzip and index resulting isec1-2 VCFs
	for i in "$iSAMPLE_DIR/01--isec1-2/"*.vcf ; do
		bgzip -c "$i" > "${i}.gz"
		tabix "${i}.gz"
	done

	# take isec of 1-2 result with 3rd file
	bcftools isec \
	-p "$iSAMPLE_DIR"/"02--isec3" \
	"$iSAMPLE_DIR/01--isec1-2/0002.vcf.gz" \
	"$FILE3"

	# rename final result from isec3
	mkdir "${iSAMPLE_DIR}/03--isec3-final" # make dir for final VCF
	BASENAME=$(basename "${FILE1%.gz}")
	SAMPLE_NAME="${BASENAME%-*}" # i.e. rootname without '-run1|2|3' suffix
	FULL_BASENAME_EXT=$(echo "${BASENAME#*.}")
	SEMI_OUT_BASENAME="${SAMPLE_NAME}.${FULL_BASENAME_EXT}"
	OUT_BASENAME="${SEMI_OUT_BASENAME%.vcf}.isec3.vcf"
	mv "${iSAMPLE_DIR}/02--isec3/0002.vcf" "${iSAMPLE_DIR}/03--isec3-final/${OUT_BASENAME}"

done
