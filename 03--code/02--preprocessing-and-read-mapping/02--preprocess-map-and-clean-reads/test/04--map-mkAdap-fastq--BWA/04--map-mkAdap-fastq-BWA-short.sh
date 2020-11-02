#!/bin/bash
#SBATCH --job-name bwaMemAlignment.sh 			### Job name
#SBATCH -o Output_Log-bwaMemAlignment.sh.log	### File to store output
#SBATCH -e Error_Log-bwaMemAlignment.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 24:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=7					### Number of cores per node requested. Max 20/node
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
# Aligns FASTQ sequence reads to a reference genome using BWA-MEM


################################################################################
# ARGUMENTS
################################################################################
# REF_GENOME_DATABASE_FOLDER_FILES_PREFIX=$1
# e.g. "/path/to/refGenomeFolder/prefix" where e.g. "prefix" would be "mm10_UCSC" corresponding to various files inside reference genome folder such as "mm10_UCSC.fa"
MKADAP_FASTQ=$1 # Interleaved FASTQ made with Picard SamToFastq
ALIGNED_BAM_HEADER=$2


################################################################################
# EXECUTION CODE
################################################################################
module load bwa

bwa mem -M -t 7 -p \
"/lustre/project/zpursell/leo/projects/project-02--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumors/02--raw-data/02--standard-reference/genome/mm10--UCSC/mm10_UCSC.fa" \
"$MKADAP_FASTQ" > "$ALIGNED_BAM_HEADER".bam


################################################################################
# END JOB
################################################################################
echo End Job
