#!/bin/bash
#SBATCH --job-name processSomatic.sh 			### Job name
#SBATCH -o Output_Log-processSomatic.sh.log	### File to store output
#SBATCH -e Error_Log-processSomatic.sh.err	### File to store error messages
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
processSomatic "$VARSCAN_SNP_FILE"

java -jar "/lustre/project/zpursell/leo/code-packages-external/VarScan.v2.4.3.jar" \
processSomatic "$VARSCAN_INDEL_FILE"


################################################################################
# END JOB
################################################################################
echo End Job
