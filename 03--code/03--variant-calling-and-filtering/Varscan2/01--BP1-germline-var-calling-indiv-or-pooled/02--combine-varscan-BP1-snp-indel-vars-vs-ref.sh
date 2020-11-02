#!/usr/bin/env zsh

# title: 02--combine-varscan-BP1-snp-indel-vars-vs-ref.sh
# author: Leo Williams
# date: 2020-07-06 7:41 PM

# ==============================================================================

# USAGE ########################################



# ABOUT ########################################

# This script takes as input a top-level folder containing folders for each sample, each of which contain a results dir with regular and filtered snp and indel variants according to the output of script '01--varscan-BP1-variants-vs-reference.sh'. It goes through each sample and combines the regular SNP and INDEL calls into an 'all' file, and does the same thing for the filtered SNP and INDEL sets.

# ARGUMENTS ########################################

IN_DIR=$1

# SETUP ########################################

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$IN_DIR"

# FUNCTIONS ########################################

	concatVCF () {
		# define variables
		GZVCF1="$1"
		GZVCF2="$2"
		OUT_VCF_PREFIX="$3"

		if [[ -z "$1" || -z "$2" || -z "$3" ]]; then
			echo "ERROR: Needs 3 arguments."
			echo "Usage: concatVCF <first.vcf.gz> <second.vcf.gz> <'output file prefix'>"
		fi

		# combine VCFs
		bcftools concat \
		--allow-overlaps \
		--rm-dups none \
		"$GZVCF1" "$GZVCF2" > "${OUT_VCF_PREFIX}.vcf"
		# --output-type v \

		# NOTE: removing .gz and .tbi files at the end anyway
		# # compress and index
		# bgzip -c "${OUT_VCF_PREFIX}.vcf" > "${OUT_VCF_PREFIX}.vcf.gz"
		# tabix "${OUT_VCF_PREFIX}.vcf.gz"
	}


# MAIN CODE ########################################

# for every VCF in the IN_DIR_OF_VARSCAN2_VCF_DIRS hierarchy, compress with bgzip and index with tabix
for i in "$IN_DIR"/**/*.vcf; do
	# echo "$i"
	bgzip -c "$i" > "${i}.gz"
	tabix "${i}.gz"
done

# for each sample--as identified by each *.snp.vcf.gz file that is present within each sample results folder:
for iSAMPLE_SNP_VCF in "$IN_DIR"/**/*.snp.vcf.gz ; do
	# echo iter name
	echo "$iSAMPLE_SNP_VCF"

	# cd to dir with VCFs for this sample
	cd "${iSAMPLE_SNP_VCF%/*}"
	echo "Current dir:"
	pwd

	# extract file rootname to use as sample name
	FILE_BASENAME=$(basename $iSAMPLE_SNP_VCF)
	echo "$FILE_BASENAME"
	ROOT_SAMPLE_NAME=${FILE_BASENAME%%.*}
	echo "$ROOT_SAMPLE_NAME"

	# use bcftools to concatenate the indel and SNP\V files into asingle *.all.vcf.gz file
	# concat all Somatic variants
	concatVCF "${ROOT_SAMPLE_NAME}.snp.vcf.gz" \
	"${ROOT_SAMPLE_NAME}.indel.vcf.gz" \
	"${ROOT_SAMPLE_NAME}.all"

	# concat all filtered variants
	concatVCF "${ROOT_SAMPLE_NAME}.snp.filter.vcf.gz" \
	"${ROOT_SAMPLE_NAME}.indel.filter.vcf.gz" \
	"${ROOT_SAMPLE_NAME}.all.filter"
done

# remove all .vcf.gz and .tbi files to leave plain VCFs
for i in "$IN_DIR"/**/*.{gz,tbi} ; do
	rm "$i"
done