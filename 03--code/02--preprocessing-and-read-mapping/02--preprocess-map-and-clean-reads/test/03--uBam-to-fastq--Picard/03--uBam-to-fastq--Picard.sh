#!/bin/bash
#SBATCH --job-name ubamToFastq.sh 			### Job name
#SBATCH -o Output_Log-ubamToFastq.sh.log	### File to store output
#SBATCH -e Error_Log-ubamToFastq.sh.err	### File to store error messages
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
# Uses Picard SamToFastq to convert sequencing adapter-marked uBAM back to FASTQ in order to assign adapters low quality so they don't get mistaken for a read

################################################################################
# ARGUMENTS
################################################################################
MKADAP_UBAM=$1
MKADAP_FASTQ_HEADER=$2


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk SamToFastq \
-INPUT="$MKADAP_UBAM" \
-FASTQ="$MKADAP_FASTQ_HEADER".fq \
-CLIPPING_ATTRIBUTE=XT \
-CLIPPING_ACTION=2 \
-INTERLEAVE=true \
-NON_PF=true \
-TMP_DIR=$TMPDIR


################################################################################
# END JOB
################################################################################
echo End Job
