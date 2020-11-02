#!/bin/bash

################################################################################
# JOB DESCRIPTION
################################################################################
# Step 3 of VarScan2 Basic Protocol 2: Somatic Mutation Calling in Tumor-Normal Pairs


################################################################################
# ARGUMENTS
################################################################################
VARSCAN_SNP_FILE=$1
VARSCAN_INDEL_FILE=$2


################################################################################
# EXECUTION CODE
################################################################################
# 3. Run the processSomatic subcommand to divide the output files into separate files based on somatic status and confidence:

java -jar "/lustre/project/zpursell/leo/code-packages-external/VarScan.v2.4.3.jar" \
processSomatic "$VARSCAN_SNP_FILE" \
--min-var-freq 0.05

java -jar "/lustre/project/zpursell/leo/code-packages-external/VarScan.v2.4.3.jar" \
processSomatic "$VARSCAN_INDEL_FILE" \
--min-var-freq 0.05


################################################################################
# END JOB
################################################################################
echo End Job
