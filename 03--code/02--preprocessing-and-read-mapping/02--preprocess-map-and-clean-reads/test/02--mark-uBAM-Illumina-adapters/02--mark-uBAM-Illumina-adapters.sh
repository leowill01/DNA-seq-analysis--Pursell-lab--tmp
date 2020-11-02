#!/bin/bash
#SBATCH --job-name markIlluminaAdapters.sh 			### Job name
#SBATCH -o Output_Log-markIlluminaAdapters.sh.log	### File to store output
#SBATCH -e Error_Log-markIlluminaAdapters.sh.err	### File to store error messages
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
# Reads a SAM or BAM file and rewrites it with new adapter-trimming tags.
# This tool clears any existing adapter-trimming tags (XT:i:) in the optional tag region of a SAM file. The SAM/BAM file must be sorted by query name.
# Outputs a metrics file histogram showing counts of bases_clipped per read.


################################################################################
# ARGUMENTS
################################################################################
IN_UBAM=$1
OUT_HEADER=$2


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk MarkIlluminaAdapters \
-I="$IN_UBAM" \
-O="$OUT_HEADER".mkadap.u.bam \
-M="$OUT_HEADER".mkadap.u.bam.metrics.txt \
-TMP_DIR=$TMPDIR

################################################################################
# END JOB
################################################################################
echo End Job
