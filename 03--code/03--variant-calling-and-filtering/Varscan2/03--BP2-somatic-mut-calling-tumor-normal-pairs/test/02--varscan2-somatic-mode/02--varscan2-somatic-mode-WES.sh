#!/bin/bash
#SBATCH --job-name somaticModeVarScan2VCF.sh 			### Job name
#SBATCH -o Output_Log-somaticModeVarScan2VCF.sh.log	### File to store output
#SBATCH -e Error_Log-somaticModeVarScan2VCF.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 24:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=5					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmark: 5 cores
# WARNING: mutect (another variant caller) timed out at 12h on 5 cores when calling on a raw BAM. But varscan seems to be significantly faster

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Step 2 of VarScan2 Basic Protocol 2: Somatic Mutation Calling in Tumor-Normal Pairs


################################################################################
# ARGUMENTS
################################################################################
MPILEUP_FILE=$1
OUTPUT_PREFIX=$2


################################################################################
# EXECUTION CODE
################################################################################
# 2. Run VarScan2 in somatic mode using the tumor-normal mpileup file and an output prefix:
java -jar "/lustre/project/zpursell/leo/code-packages-external/VarScan.v2.4.3.jar" \
somatic "$MPILEUP_FILE" \
"$OUTPUT_PREFIX" \
--mpileup 1 \
--min-coverage 10 \
--min-var-freq 0.05 \
--somatic-p-value 0.05 \
--output-vcf 1

# TODO: --strand-filter 1 ?


################################################################################
# END JOB
################################################################################
echo End Job
