#!/bin/bash

################################################################################
# JOB DESCRIPTION
################################################################################
# Step 4 of VarScan2 Basic Protocol 2: Somatic Mutation Calling in Tumor-Normal Pairs


################################################################################
# ARGUMENTS
################################################################################
VARSCAN_SNP_HC_FILE=$1
VARSCAN_INDEL_FILE=$2
OUT_FILE_PREFIX=$3

################################################################################
# EXECUTION CODE
################################################################################
# 4. Use the somaticFilter subcommand to identify and remove somatic mutations that are likely false positives due to alignment problems near indels.
# After this step, the protocol says to run Support Protocol 1 to further filter false positives, but one of the tools for that only works on linux systems and not macOS.

# NOTE: For this script, make sure input is VCF

java -jar "/lustre/project/zpursell/leo/code-packages-external/VarScan.v2.4.3.jar" \
somaticFilter "$VARSCAN_SNP_HC_FILE" \
--indel-file "$VARSCAN_INDEL_FILE" \
--output-file "$OUT_FILE_PREFIX".somaticFilter.vcf \
--min-var-freq 0.05


################################################################################
# END JOB
################################################################################
echo End Job
