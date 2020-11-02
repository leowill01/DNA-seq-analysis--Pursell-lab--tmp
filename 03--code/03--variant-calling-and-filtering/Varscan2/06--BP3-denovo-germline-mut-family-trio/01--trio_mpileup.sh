#!/bin/bash
#SBATCH --job-name trio.sh 			### Job name
#SBATCH -o Output_Log-trio.sh.log	### File to store output
#SBATCH -e Error_Log-trio.sh.err	### File to store error messages
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
# Use VarScan2 to call de novo germline variants on family trio samples

################################################################################
# ARGUMENTS
################################################################################
REF_FASTA=$1
FATHER_BAM=$2
MOTHER_BAM=$3
CHILD_BAM=$4
TRIO_HEADER=$5


################################################################################
# EXECUTION CODE
################################################################################
module load samtools

samtools mpileup -B -q 1 \
-f $REF_FASTA \
$FATHER_BAM \
$MOTHER_BAM \
$CHILD_BAM \
> "$TRIO_HEADER".mpileup

################################################################################
# END JOB
################################################################################
echo End Job
