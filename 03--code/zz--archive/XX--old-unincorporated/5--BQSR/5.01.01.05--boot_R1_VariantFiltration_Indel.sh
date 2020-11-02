#!/bin/bash
#SBATCH --job-name variantFiltrationIndel.sh 			### Job name
#SBATCH -o Output_Log-variantFiltrationIndel.sh.log		### File to store output
#SBATCH -e Error_Log-variantFiltrationIndel.sh.err		### File to store error messages
#SBATCH --qos normal								### Quality of service queue
#SBATCH -t 12:00:00									### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1									### Nodes requested for job
#SBATCH --ntasks-per-node=1							### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1							### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G									### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmark: 1 core
# NOTE: Run interactively to avoid queue time

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Filter high-confidence variants from a raw bootstrapped germIndel.vcf generated with HaplotypeCaller to make a knownSites.vcf training set to feed back into the BaseRecalibrator algorithm. This process will be repeated until convergence is obtained (AnalyzeCovariates).

################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1
IN_RAW_INDEL_VCF=$2
OUT_HARD_FILTERED_INDEL_VCF_HEADER=$3
# -O="HEADER".hardFilter.vcf


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk VariantFiltration \
-R="$REF_GENOME_FASTA" \
-V="$IN_RAW_INDEL_VCF" \
-O="$OUT_HARD_FILTERED_INDEL_VCF_HEADER".hardFilter.vcf \
--filter-expression="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0  || SOR > 10.0" \
--filter-name="IndelHardFilter" \
-TMP_DIR=$TMPDIR

################################################################################
# END JOB
################################################################################
echo End Job
