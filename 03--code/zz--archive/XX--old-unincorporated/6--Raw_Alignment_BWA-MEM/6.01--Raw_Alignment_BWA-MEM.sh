#!/bin/bash
#SBATCH --job-name rawAlignBwaMem.sh 			### Job name
#SBATCH -o Output_Log-rawAlignBwaMem.sh.log	### File to store output
#SBATCH -e Error_Log-rawAlignBwaMem.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=10					### Number of cores per node requested. Max 20/node
#SBATCH --mem=10G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmarked for 10 cores, 10GB memory

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Aligns a sample's fwd and rev FASTQ read files to a reference genome FASTA
# The tool is run similar to the options when realigning the interleaved FASTQ during preprocessing except without the -p flag since the inputs are F/R FASTQ files and not a single interleaved FASTQ.


################################################################################
# ARGUMENTS
################################################################################
REF_GENOME_FASTA=$1
FWD_READS_FASTQ=$2
REV_READS_FASTQ=$3
OUT_BAM_PREFIX=$4


################################################################################
# EXECUTION CODE
################################################################################
module load bwa

bwa mem -M -t 10 \
"$REF_GENOME_FASTA" \
"$FWD_READS_FASTQ" \
"$REV_READS_FASTQ" \
> "$OUT_BAM_PREFIX".raw.bam

# The FASTA filepath is actually the prefix for the index files e.g. "mm10_UCSC.fa.amb"


################################################################################
# END JOB
################################################################################
echo End Job
