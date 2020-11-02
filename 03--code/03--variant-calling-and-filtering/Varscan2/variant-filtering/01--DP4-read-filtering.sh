#!/usr/bin/env zsh

# title: filter-only-tumor-reads.sh
# author: Leo Williams
# date: 2020 08_Aug 02

# USAGE ########################################

#

# ABOUT ########################################

# This script takes a recursive dir with VCFs (tumor-normal somatic mutation VCFs), finds every VCF, and uses bcftools to filter only for variants that have 1. no reads in the normal sample, and 2. at least 1 read in each direction in the tumor sample. It outputs the filtered callset to a new VCF.

# WARNING: The VCF files used to write this script are output from VarScan2 Basic Protocol 2 tumor-normal somatic mutation calling. The line of text for each variant ends with the DP4 values for the tumor, which code in this script is specific to. This may not work with VCF outputs from other programs.

# ARGUMENTS ########################################

IN_DIR_VCFS=$1
OUT_DIR=$2

# SETUP ########################################

# Make timestamped dir for results output
RESULTS_DIR="${OUT_DIR}/results--filterReads-$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir $RESULTS_DIR

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$RESULTS_DIR"

# Make tmp dir to hold intermediate files
DIR_TMP="$RESULTS_DIR"/tmp
mkdir "$DIR_TMP"

# MAIN CODE ########################################

# for every VCF in the IN_DIR_OF_VCFS hierarchy, compress with bgzip and index with tabix
for i in "$IN_DIR_VCFS"/**/*.vcf; do
	echo "$i" # print filename
	iBASENAME=$(basename "$i") # get basename & remove path
	bgzip -c "$i" > "${DIR_TMP}/${iBASENAME}.gz"
	tabix "${DIR_TMP}/${iBASENAME}.gz"
done

for gzvcf in "$DIR_TMP"/*.vcf.gz ; do
	# get file basename for output file
	gzBASENAME=$(basename "$gzvcf")
	echo $gzBASENAME

	# remove variants that with: 1) any alt-reads in the normal, 2) at least 1 read in both directions in the tumor
		# define output filename
		gzvcfFiltDP4="${gzBASENAME%.vcf.gz}.filtDP4.vcf"

		# filter variants
		bcftools view \
		--include 'AD[0] == 0 && DP4[1:2] > 0 && DP4[1:3] > 0' \
		"$gzvcf" \
		> "${RESULTS_DIR}/${gzvcfFiltDP4}"

		# # bgzip and index normal read-filtered VCF
		# bgzipAndTabixVcf "$tmpVcfFiltNReads"

	# TODO: commenting out to implement suggestions at [https://github.com/samtools/bcftools/issues/1275#issuecomment-674798948]
	# # take the normal reads-filtered .vcf.gz file and filter that for only variants that have at least 1 read in both directions in the tumor. An inverse grep is used because bcftools cannot filter based on sample-specific DP4 values.
	# 	# define output filename
	# 	outVcf_FiltNReads_FiltTReads="${RESULTS_DIR}/${tmpVcfFiltNReads%.filtNReads.vcf}.filtDP4.vcf"

	# 	# filter variants
	# 	grep -v \
	# 	-e '0,.$' \
	# 	-e '.,0$' \
	# 	"$DIR_TMP/$tmpVcfFiltNReads" \
	# 	> "$outVcf_FiltNReads_FiltTReads"
done

rm -rf "$DIR_TMP"

# NOTE: implemented with function countVcfVariants() in .zshrc
# # print a count of number of variants in each sample
# 	# make counts summary file
# 	COUNTS_FILE="${RESULTS_DIR}/variant-counts.txt"
# 	printf "" > "$COUNTS_FILE"

# 	# add counts for each sample
# 	for vcf in "$RESULTS_DIR"/*.vcf ; do
# 		echo $(basename "$vcf") >> "$COUNTS_FILE"
# 		grep '^[^#]' "$vcf" | wc -l >> "$COUNTS_FILE"
# 	done

# ===============================

	# bcftools view \
	# --include 'FORMAT/AD[0] == 0 & DP4[1:2] > 0' \
	# "$gz" > "$RESULTS_DIR"/"${gzBASENAME%.vcf.gz}".readFilter.vcf
	# # --include 'FORMAT/AD[0] == 0 & DP4[1:2] > 0 & DP4[1:3] > 0' \
	# # --include 'FORMAT/AD[0] == 0 && FORMAT/AD[1] >= 2' \
	# # --exclude 'DP4[1:2] == 0 | DP4[1:3] == 0' \

	# bcftools view \
	# --include 'FORMAT/AD[0] == 0' \
	# "$gz" \
	# | \
	# bcftools view \
	# --exclude 'DP4[1:2] == 0' \
	# - > "$RESULTS_DIR"/"${gzBASENAME%.vcf.gz}".readFilter.vcf
	# # --include 'FORMAT/AD[0] == 0 && FORMAT/AD[1] >= 2' \
	# # --exclude 'DP4[1:2] == 0 | DP4[1:3] == 0' \


	# --include 'DP4[1:2] >= 1' \
	# --include 'DP4[1:2] >= 1 && DP4[1:3] >= 1' \
	# --exclude 'DP4[2] == 0 && '

# rm -rf "$DIR_TMP"