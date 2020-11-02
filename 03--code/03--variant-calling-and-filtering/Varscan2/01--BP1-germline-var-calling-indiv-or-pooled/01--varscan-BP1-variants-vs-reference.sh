#!/usr/bin/env zsh

# FILENAME:		varscan-germline-variants-vs-reference.sh
# AUTHOR:		Leo Williams
# DATE:			2020-04-28
# USAGE:		[USAGE]
# DESCRIPTION:	This script calls germline variant (SNP) differences in a single vs. a reference genome. Code is adapted from the VarScan2 Basic Protocol 1: Germline Variant Calling in Individual or Pooled Samples [Koboldt et al., 2013]
# VERSION:		[VERSION]
# NOTES:		[NOTES]
# ZSH VERSION:	[ZSH VERSION]
# DEV PLATFORM:	macOS 10.15.4

# ==============================================================================

# ARGUMENTS ########################################

REF_FASTA=$1
IN_BAM=$2
OUTPUT_DIR=$3

# SETUP ########################################

# Make timestamped dir for results output
RESULTS_DIR="${OUTPUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir $RESULTS_DIR

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$RESULTS_DIR"

# MAIN CODE ########################################

# Basic Protocol 1: Germline Variant Calling in Individual or Pooled Samples
# Detection of SNVs and indels in an individual sample is a common task for NGS studies. With VarScan, it is possible to call both types of variants simultaneously using the mpileup output from SAMtools. At each position, VarScan looks for variants that meet user-defined minimum criteria for sequencing coverage, number of supporting reads, variant allele frequency (VAF), and Fisher's Exact Test p-value.
# The subcommands for germline variant calling are mpileup2snp (for SNV calling), mpileup2indel (for indel calling), and mpileup2cns (for consensus or simultaneous SNV/indel calling).
# Relevant user-defined parameters are listed in Table 2. Among these, the most important are the --min-coverage and --min-var-freq parameters, which essentially govern VarScan's sensitivity and specificity for variant calling. Conservative values ( --mincoverage 20 --min-var-freq 0.20) yield high-confidence calls, but may miss variants that are under-covered or under-sampled in the dataset. For pooled samples, one should generally specify higher minimum coverage but lower variant allele frequency thresholds.

# call variants from a BAM file

	# 1. run samtools mpileup on the BAM file
	OUT_MPILEUP="${RESULTS_DIR}/${OUT_PREFIX}.mpileup"
	samtools mpileup -B -q 1 -f "$REF_FASTA" "$IN_BAM" > "$OUT_MPILEUP"

	# 2. run varscan mpileup2snp to call SNPs
	OUT_SNP_VCF="${RESULTS_DIR}/${OUT_PREFIX}.snp.vcf"
	varscan mpileup2snp \
	"$OUT_MPILEUP" \
	--min-coverage 10 \
	--min-var-freq 0.20 \
	--p-value 0.05 \
	--output-vcf 1 \
	> "$OUT_SNP_VCF"

	# 3. run varscan mpileup2indel to call indels
	OUT_INDEL_VCF="${RESULTS_DIR}/${OUT_PREFIX}.indel.vcf"
	varscan mpileup2indel \
	"$OUT_MPILEUP" \
	--min-coverage 10 \
	--min-var-freq 0.10 \
	--p-value 0.10 \
	--output-vcf 1 \
	> "$OUT_INDEL_VCF"
	# Note that we use slightly less stringent parameters for indel calling, which reflects the fact that these variants are more difficult to detect in NGS data. VarScan's sensitivity and range for indel detection are dictated by the underlying read alignments. For example, aligning paired-end 100-bp reads with the BWA algorithm typically allows for detection of indels in the 1-30 bp size range, corresponding to the gap size that the aligner allows. Alternate indel detection strategies, such as the split-read mapping approach employed by Pindel (Ye, Schulz et al. 2009), are recommended for larger indels.

		# remove mpileup file
		rm "$OUT_MPILEUP"

	# 4. filter SNP calls to remove those near indel positions
	OUT_SNP_FILTER_VCF="${RESULTS_DIR}/${OUT_PREFIX}.snp.filter.vcf"
	varscan filter \
	"$OUT_SNP_VCF" \
	--indel-file "$OUT_INDEL_VCF" \
	--output-file "$OUT_SNP_FILTER_VCF" \
	--output-vcf 1

	# 5. filter indels to obtain a higher confidence set
	OUT_INDEL_FILTER_VCF="${RESULTS_DIR}/${OUT_PREFIX}.indel.filter.vcf"
	varscan filter \
	"$OUT_INDEL_VCF" \
	--min-reads2 4 \
	--min-var-freq 0.15 \
	--p-value 0.05 \
	--output-vcf 1 \
	--output-file "$OUT_INDEL_FILTER_VCF"

	# END protocol instructions

# # Filter for variants that are germline as determined by VAF, i.e. VAF=[0.40-0.60]
# varscan filter \
# "${RESULTS_DIR}/${OUT_PREFIX}.snp.filter.vcf" \
# --min-var-freq 0.40 \
# --min-strands2 2 \
# --indel-file "${RESULTS_DIR}/${OUT_PREFIX}.indel.vcf" \
# --output-file "${RESULTS_DIR}/${OUT_PREFIX}.snp.filter.minVaf40.vcf"

# varscan filter \
# "${RESULTS_DIR}/${OUT_PREFIX}.indel.filter.vcf" \
# --min-var-freq 0.40 \
# --min-strands2 2 \
# --indel-file "${RESULTS_DIR}/${OUT_PREFIX}.indel.vcf" \
# --output-file "${RESULTS_DIR}/${OUT_PREFIX}.snp.filter.minVaf40.vcf"


# See Support Protocol 1 for additional filtering recommendations