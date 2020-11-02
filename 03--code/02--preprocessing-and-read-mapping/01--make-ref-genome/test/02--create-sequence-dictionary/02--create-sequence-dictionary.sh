#!/bin/bash
#SBATCH --job-name refGenomeDictPicard.sh 			### Job name
#SBATCH -o Output_Log-refGenomeDictPicard.sh.log	### File to store output
#SBATCH -e Error_Log-refGenomeDictPicard.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=40G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmarked: 1 CPU core, 40GB memory (25s runtime)

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Creates a dictionary for the reference genome for downstream use of other Picard & GATK4 tools.

################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1


################################################################################
# EXECUTION CODE
################################################################################
# Generate reference genome sequence dictionary with Picard (GATK)
module load gatk
gatk CreateSequenceDictionary \
-R "$REF_GENOME_FASTA" \
-TMP_DIR=$TMPDIR


################################################################################
# END JOB
################################################################################
echo End Job
