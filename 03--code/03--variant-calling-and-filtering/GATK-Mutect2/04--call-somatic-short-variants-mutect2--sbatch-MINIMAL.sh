#!/usr/bin/env bash

#SBATCH --job-name=mutect2TN    					### Job name
#SBATCH --output=Output_Log-mutect2TN.log		### File to store output
#SBATCH --error=Error_Log-mutect2TN.err			### File to store error messages
#SBATCH --qos=normal				    	### Quality of service queue
#SBATCH --time=24:00:00						### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1 				### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwilli24@tulane.edu

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

# setup
module load gatk/4.1.8.1

# args
in_ref_fasta=$1
in_tumor_bam=$2
in_normal_bam=$3
normal_bam_rg_sample_name=$4
out_vcf_gz=$5

# run mutect
gatk Mutect2 \
--reference "$in_ref_fasta" \
--input "$in_tumor_bam" \
--input "$in_normal_bam" \
--normal-sample "$normal_bam_rg_sample_name" \
--output "$out_vcf_gz"

# filter the raw Mutect2 calls
gatk FilterMutectCalls \
--reference "$in_ref_fasta" \
--variant "$out_vcf_gz" \
--output "$out_vcf_gz_filtered"
# --contamination-table contamination.table \
# --tumor-segmentation segments.tsv \

# END ########################################
echo "End"