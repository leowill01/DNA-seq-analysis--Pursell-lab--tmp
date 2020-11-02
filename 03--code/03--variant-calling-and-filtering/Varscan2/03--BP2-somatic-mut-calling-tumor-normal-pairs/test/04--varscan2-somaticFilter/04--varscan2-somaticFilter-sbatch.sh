#!/bin/bash
#SBATCH --job-name somaticFilter.sh 			### Job name
#SBATCH -o Output_Log-somaticFilter.sh.log	### File to store output
#SBATCH -e Error_Log-somaticFilter.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 24:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE benchmarked
# NOTE Run interactively ro reduce queue time


echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

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
--output-file "$OUT_FILE_PREFIX".somaticFilter.vcf


################################################################################
# END JOB
################################################################################
echo End Job
