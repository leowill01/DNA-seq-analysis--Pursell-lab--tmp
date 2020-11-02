#!/usr/bin/env zsh

# title:		02--Varscan-AP1-combineCohortGermSNV_InDel.sh
# author:		Leo Williams
# date:			2020-07-08 4:31 PM

# USAGE ########################################

	# zsh combine-Varscan2-variants.sh $IN_DIR_OF_VARSCAN_AP1_VCFS

	# Run from within script dir and put full filepath of input dir as arg

# ABOUT ########################################

	# This script searches for results from VarScan2 Alternate Protocol 1: Germline Variant Calling in a Cohort of Individuals and combines the 'indel' and 'snp' files into an additional '*.all.vcf' VCF file containing all indels and SNVs for the cohort.

	# This script works on files in-place and doesn't require a results output directory.

# ARGUMENTS ########################################

	IN_DIR_OF_VARSCAN_AP1_VCFS=$1 # can be recursive with multiple sample folders

# SETUP ########################################

	# Make a copy of the script run and put it into the results dir
	SCRIPT_PATH=$(echo "$0")
	cp "$SCRIPT_PATH" "$IN_DIR_OF_VARSCAN_AP1_VCFS"

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

	# for every VCF in the IN_DIR_OF_VARSCAN_AP1_VCFS hierarchy, compress with bgzip and index with tabix
	for i in "$IN_DIR_OF_VARSCAN_AP1_VCFS"/**/*.vcf; do
		# echo "$i"
		bgzip -c "$i" > "${i}.gz"
		tabix "${i}.gz"
	done

	# for each sample--as identified by each *.snp.vcf.gz file that is present within each sample results folder:
	for iSAMPLE_SNP_VCF in "$IN_DIR_OF_VARSCAN_AP1_VCFS"/**/*.snp.vcf.gz ; do
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

		# use bcftools to concatenate the indel and SNP\V files respectively for each {Somatic|Germline|LOH} class into single *.{Somatic|Germline|LOH}.all.vcf.gz files
		# concat all cohort SNP and InDel variants
		concatVCF "${ROOT_SAMPLE_NAME}.cohort.varScan.snp.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.cohort.varScan.indel.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.cohort.varScan.all"
	done

	# remove all .vcf.gz and .tbi files to leave plain VCFs
	for i in "$IN_DIR_OF_VARSCAN_AP1_VCFS"/**/*.{gz,tbi} ; do
		rm "$i"
	done
