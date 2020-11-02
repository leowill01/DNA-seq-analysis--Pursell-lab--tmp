#!/usr/bin/env zsh

# FILENAME:		04--filter-somatic-vars-high-low-VAF.sh
# AUTHOR:		Leo Williams
# DATE:			Created: 2020 04_Apr 03
# USAGE:		[USAGE]
# DESCRIPTION:	See https://stackoverflow.com/questions/430078/shell-script-templates for template ideas
# VERSION:		[VERSION]
# NOTES:		[NOTES]
# ZSH VERSION:	[ZSH VERSION]
# DEV PLATFORM:	[DEV PLATFORM]

# ==============================================================================

# ARGUMENTS ########################################

# INPUT=$1
# OUTPUT_DIR=$2
INPUT="/Volumes/Pursell-Lab-HD/Pursell-lab/projects/project_02--DNA-seq-analysis/04--analysis/06--mutation-analysis/09--HCT116-Selina-revisions/02--somatic-variants/00--input-vcfs"
OUTPUT_DIR="/Volumes/Pursell-Lab-HD/Pursell-lab/projects/project_02--DNA-seq-analysis/04--analysis/06--mutation-analysis/09--HCT116-Selina-revisions/02--somatic-variants/04--filter-by-VAF"

# SETUP ########################################

# Make timestamped dir for results output
RESULTS="${OUTPUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir $RESULTS

# MAIN CODE ########################################

# pseudocode
# for every **/*.vcf in dir, 
# 	- bgzip and tabix (then remove when done)
# 	- bcftools view --include ''

for i in "$INPUT"/**/*.vcf; do
	bgzip -c "$i" > "${i}.gz"
	tabix "${i}.gz"
done

for i in "$INPUT"/**/*.vcf.gz; do
	bcftools view \
	--include 'FORMAT/FREQ[1] > '
done

# TEST ########################################

# SETUP
test_vcf="S459F-6K2_T-6K_N.all.Somatic.hc.filter.hg38_multianno.vcf"
test_gz="S459F-6K2_T-6K_N.all.Somatic.hc.filter.hg38_multianno.vcf.gz"
# test
bgzip -c $test_vcf > "${test_vcf}.gz"
tabix "${test_vcf}.gz"

# try to filter by varscan FREQ
bcftools view --include 'FORMAT/FREQ[1] > 40.00 & FORMAT/FREQ[1] < 60.00' $test_gz > 01-test-out.vcf # FIXME: DOES NOT WORK

# extract all vars that have any DP4 value of 0
bcftools query --print-header --include 'DP4 = 0' -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP4]\n' $test_gz > 02-test-out.vcf
# extract all vars that have any DP4 value of 0 in the normal
bcftools query --print-header --include 'DP4[0:*] = 0' -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP4]\n' $test_gz > 02-test-out.vcf # BUG: doesnt work - extracts normal vars where no DP4 values are 0
# extract any vars with DP4 4th value = 0
bcftools query --print-header --include 'DP4[*:3] = 0' -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP4]\n' $test_gz > 02-test-out.vcf # BUG: still extracts vars where DP4 4th value is not 0

# extract vars with on-the-gly calc AF > 0.5
bcftools query --include 'AF >= 0.5' -f '%CHROM\t%POS\t%REF\t%ALT[\t%FREQ\t%DP4]\n' $test_gz > 03-test-out.vcf # BUG: doesn't work bc manual filtering shows many more variants

# try with manual VAF DP4 calculation for tumor sample (sample '1' on 0-index)
bcftools query --print-header --include '((DP4[2:2] + DP4[2:3]) / (SUM(DP4[2]))) >= 0.5' -f '%CHROM\t%POS\t%REF\t%ALT[\t%FREQ\t%DP4]\n' $test_gz > 04-test-out.vcf


# bcftools query --print-header --include '(DP4[1:2]+DP4[1:3])/(sum(DP4[1])) >= 0.4' -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP4]\n' $test_gz > 02-test-out.vcf
# bcftools view --include 'FORMAT/FREQ > 40' $test_gz > 03-test-out.vcf
# bcftools view --include 'FORMAT/DP4[1:3] = 0' $test_gz > 03-test-out.vcf
# bcftools view --include '(DP4[1:2]+DP4[1:3])/(DP4[1:0]+DP4[1:1]+DP4[1:2]+DP4[1:3]) >= 0.4' $test_gz > 03-test-out.vcf
# bcftools view -i'(DP4[0:0]+DP4[0:1])/(DP4[0:2]+DP4[0:3]) > 0.3' $test_gz > 03-test-out.vcf
# bcftools view --include '((FORMAT/DP[1]) > 40.00) &  < 60.00' $test_gz > 04-test-out.vcf
# bcftools view --include '((FORMAT/AD[1])/(FORMAT/DP[1])) > 0.40 & ((FORMAT/AD[1])/(FORMAT/DP[1])) < 0.60' $test_gz > 03-test-out.vcf
# query just FREQ to see if they match up
# bcftools query --print-header --include '((FORMAT/AD[2])/(FORMAT/DP[2])) > 0.40 & ((FORMAT/AD[2])/(FORMAT/DP[2])) < 0.60' -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD\t%FREQ]\n' $test_gz > 04-test-out.vcf
# bcftools query --print-header --include '((FORMAT/AD[2])/(FORMAT/DP[2])) > 0.40 & ((FORMAT/AD[2])/(FORMAT/DP[2])) < 0.60' -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD\t%FREQ]\n' $test_gz > 04-test-out.vcf


# VCFTOOLS ########################################

# vcftools [ --vcf FILE | --gzvcf FILE | --bcf FILE] [ --out OUTPUT PREFIX ] [ FILTERING OPTIONS ] [ OUTPUT OPTIONS ]

# manual examples

	# Output allele frequency for all sites in the input vcf file from chromosome 1
	vcftools --vcf $test_vcf --freq --chr 1 --out chr1_analysis