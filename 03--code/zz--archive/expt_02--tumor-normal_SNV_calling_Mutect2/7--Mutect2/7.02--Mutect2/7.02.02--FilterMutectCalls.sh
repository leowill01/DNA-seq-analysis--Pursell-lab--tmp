#!/bin/bash
#SBATCH --job-name filterMutectCalls.sh 			### Job name
#SBATCH -o Output_Log-filterMutectCalls.sh.log	### File to store output
#SBATCH -e Error_Log-filterMutectCalls.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmark: 1 core, 64GB memory (10s runtime)

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Takes raw Mutect2 somatic mutation output and filters only the high-confidence calls from that VCF dataset

################################################################################
# ARGUMENTS
################################################################################
IN_RAW_MUTECT_VARIANTS_VCF=$1
OUT_MUTECT_HC_SOMATIC_MUTATIONS_VCF_HEADER=$2

################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk FilterMutectCalls \
-V "$IN_RAW_MUTECT_VARIANTS_VCF" \
-O "$OUT_MUTECT_HC_SOMATIC_MUTATIONS_VCF_HEADER".mutectFilter.vcf \

################################################################################
# END JOB
################################################################################
echo End Job
