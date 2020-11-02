#!/bin/bash
#SBATCH --job-name MergeBamAlignment.sh 			### Job name
#SBATCH -o Output_Log-MergeBamAlignment.sh.log	### File to store output
#SBATCH -e Error_Log-MergeBamAlignment.sh.err	### File to store error messages
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
# Use GATK MergeBamAlignment to add metadata from uBAM back to the adapter-marked/removed/discounted and aligned BAM and to sort the final BAM by coordinate, which is required for use by MarkDuplicates

################################################################################
# ARGUMENTS
################################################################################
MKADAP_UBAM=$1 # from FastqToSam
MKADAP_BAM=$2 # from BWA-MEM
OUT_MERGED_BAM_HEADER=$3


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk MergeBamAlignment \
-R="/lustre/project/zpursell/leo/projects/project-02--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumors/02--raw-data/02--standard-reference/genome/mm10--UCSC/mm10_UCSC.fa" \
-UNMAPPED="$MKADAP_UBAM" \
-ALIGNED="$MKADAP_BAM" \
-O="$OUT_MERGED_BAM_HEADER".merged.sorted.bam \
--SORT_ORDER=coordinate \
--CREATE_INDEX=true \
--ADD_MATE_CIGAR=true \
--CLIP_ADAPTERS=false \
--CLIP_OVERLAPPING_READS=true \
--INCLUDE_SECONDARY_ALIGNMENTS=true \
--MAX_INSERTIONS_OR_DELETIONS=-1 \
--PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
--ATTRIBUTES_TO_RETAIN=XS \
-TMP_DIR=$TMPDIR



################################################################################
# END JOB
################################################################################
echo End Job
