#!/usr/bin/env zsh

# title: 03--maskAP1CohortGermVarsFromTumor_FauxnormalSomaticMuts.sh
# author: Leo Williams
# date: 2020-07-08 7:43 PM

# USAGE ########################################

# zsh 03--maskAP1CohortGermVarsFromTumor_FauxnormalSomaticMuts.sh \
# $IN_DIR_OF_VCFS \
# $IN_COHORT_GERM_VARS_VCF_GZ \
# $OUT_DIR

# ABOUT ########################################

# This script uses output from Varscan's Alternate Protocol 1: Germline Variant Calling in a Cohort of Individuals and Basic Protocol 2: Somatic Mutation Detection in Tumor-Normal Pairs. This script is used when no matched normal was available for a cohort of related samples and instead somatic mutations were called in BP2 using other tumors as the "normals". This results in a lot of false positives where normal bases are called as mutations

# ARGUMENTS ########################################

	IN_DIR_OF_VCFS=$1
	IN_COHORT_GERM_VARS_VCF_GZ=$2 # compressed with bgzip and indexed with tabix
	OUT_DIR=$3

# SETUP ########################################

	# Make timestamped dir for results output
	RESULTS_DIR="${OUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
	mkdir "$RESULTS_DIR"

	# Make a copy of the script run and put it into the results dir
	SCRIPT_PATH=$(echo "$0")
	cp "$SCRIPT_PATH" "$RESULTS_DIR"

# MAIN CODE ########################################

for iVCF in "$IN_DIR_OF_VCFS"/*.vcf ; do
	# print iter name
	printf "\n\n"
	printf "Working on file: ${iVCF}"

	# get VCF basename and rootname
	iVCF_BASENAME=$(basename "$iVCF")
	printf "\nBasename: ${iVCF_BASENAME}"
	iVCF_ROOTNAME=${iVCF_BASENAME%%.*}
	printf "\nRootname: ${iVCF_ROOTNAME}"

	# make sample out dir in results dir
	iVCF_OUT_DIR="${RESULTS_DIR}/${iVCF_ROOTNAME}"
	printf "\nVCF output dir: ${iVCF_OUT_DIR}"
	mkdir "$iVCF_OUT_DIR"

	# compress VCF with bgzip and index with tabix
	iVCF_GZ="$iVCF_OUT_DIR"/"${iVCF_BASENAME}.gz"
	printf "\nVCF .gz file: ${iVCF_GZ}"
	bgzip -c "$iVCF" > "$iVCF_GZ"
	tabix "$iVCF_GZ"

	# get the intersection and unique variants between the sample somatic VCF (with lots of false positives) and the cohort germline variant VCF
	ISEC_OUT_DIR="${iVCF_OUT_DIR}/isec-out"
	bcftools isec \
	-p "$ISEC_OUT_DIR" \
	"$iVCF_GZ" \
	"$IN_COHORT_GERM_VARS_VCF_GZ"

	# copy & rename VCF of unique mutations to somatic mutation VCF where cohort variants have been masked out to the sample out dir one level higher
	cp "$ISEC_OUT_DIR/0000.vcf" \
	"${iVCF_OUT_DIR}/${iVCF_BASENAME%.vcf}.masked.vcf"
done

# unzip the somatic VCFs used for easy filesize comparison with the size of the uniques file
for iGZ in "$RESULTS_DIR"/**/*.gz ; do
	gunzip -c "$iGZ" > "${iGZ%.gz}"
done

# remove all .gz and .tbi files from results
rm "$RESULTS_DIR"/**/*.{gz,tbi}