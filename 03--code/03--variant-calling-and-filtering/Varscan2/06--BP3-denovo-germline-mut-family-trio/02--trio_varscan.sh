#!/bin/bash
#SBATCH --job-name trio_02.sh 			### Job name
#SBATCH -o Output_Log-trio_02.sh.log	### File to store output
#SBATCH -e Error_Log-trio_02.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 24:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Varscan2 basic protocol 3: germline variant calling of family trios
# step 2: mpileup2snp


################################################################################
# ARGUMENTS
################################################################################
TRIO_MPILEUP=$1
TRIO_OUTPUT_PREFIX=$2


################################################################################
# EXECUTION CODE
################################################################################
java -jar /lustre/project/zpursell/leo/code_packages/VarScan.v2.4.3.jar trio \
$TRIO_MPILEUP \
"$TRIO_OUTPUT_PREFIX" \
--min-coverage 10 \
--min-var-freq 0.20 \
--p-value 0.05 \
--adj-var-freq 0.05 \
--adj-p-value 0.15


################################################################################
# END JOB
################################################################################
echo End Job
