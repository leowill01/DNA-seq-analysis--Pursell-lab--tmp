#!/usr/bin/env zsh

# FILENAME:		combine-Varscan2-variants.sh
# AUTHOR:		Leo Williams
# DATE:			2020-06-14
# VERSION:		v0.1
# ZSH VERSION:	zsh v5.8
# DEV PLATFORM:	macOS Catalina v10.15.5

# ==============================================================================

# USAGE ########################################

	# zsh combine-Varscan2-variants.sh $IN_RECUR_DIR_OF_VARSCAN_VCFS

	# Run from within script dir and put full filepath of input dir as arg

# ABOUT ########################################

	# This script searches for results from VarScan2 Basic Protocol 2: Somatic Mutation Detection in Tumor-Normal Pairs and combines the 'indel' and 'snp' files each for 'Somatic', 'LOH', and 'Germline', respectively, into additional e.g. '*.{Somatic|LOH|Germline}.all.vcf' VCF files containing all indels and SNVs for each category.

	# This script works on files in-place and doesn't require a results output directory.

# ARGUMENTS ########################################

	IN_DIR_OF_VARSCAN2_VCF_DIRS=$1 # can be recursive with multiple sample folders

# SETUP ########################################

	# Make a copy of the script run and put it into the results dir
	SCRIPT_PATH=$(echo "$0")
	cp "$SCRIPT_PATH" "$RESULTS_DIR"

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
	for i in "$IN_DIR_OF_VARSCAN2_VCF_DIRS"/**/*.vcf; do
		# echo "$i"
		bgzip -c "$i" > "${i}.gz"
		tabix "${i}.gz"
	done

	# for each sample--as identified by each *.snp.vcf.gz file that is present within each sample results folder:
	for iSAMPLE_SNP_VCF in "$IN_DIR_OF_VARSCAN2_VCF_DIRS"/**/*.snp.vcf.gz ; do
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
		# concat all Somatic variants
		concatVCF "${ROOT_SAMPLE_NAME}.snp.Somatic.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.indel.Somatic.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.all.Somatic"

		# concat all Somatic hc variants
		concatVCF "${ROOT_SAMPLE_NAME}.snp.Somatic.hc.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.indel.Somatic.hc.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.all.Somatic.hc"

		# concat all Somatic hc variants with indel-filtered hc SNV variants
		concatVCF "${ROOT_SAMPLE_NAME}.snp.Somatic.hc.filter.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.indel.Somatic.hc.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.all.Somatic.hc"

		# concat all Germline variants
		concatVCF "${ROOT_SAMPLE_NAME}.snp.Germline.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.indel.Germline.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.all.Germline"

		# concat all Germline high-confidence variants
		concatVCF "${ROOT_SAMPLE_NAME}.snp.Germline.hc.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.indel.Germline.hc.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.all.Germline.hc"

		# concat all LOH variants
		concatVCF "${ROOT_SAMPLE_NAME}.snp.LOH.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.indel.LOH.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.all.LOH"

		# concat all LOH high-confidence variants
		concatVCF "${ROOT_SAMPLE_NAME}.snp.LOH.hc.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.indel.LOH.hc.vcf.gz" \
		"${ROOT_SAMPLE_NAME}.all.LOH.hc"
	done

	# remove all .vcf.gz and .tbi files to leave plain VCFs
	for i in "$IN_DIR_OF_VARSCAN2_VCF_DIRS"/**/*.{gz,tbi} ; do
		rm "$i"
	done
