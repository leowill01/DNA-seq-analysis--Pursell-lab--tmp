#!/usr/bin/env bash

# title: Mutect2 tumor-only mode for PON creation
# author: Leo Williams
# date: 2020-09-18-17:18:54--Friday, Sep 18, 2020 05:18 PM

#SBATCH --job-name=mutect2TumorOnlyMode    					### Job name
#SBATCH --output=Output_Log-mutect2TumorOnlyMode.log		### File to store output
#SBATCH --error=Error_Log-mutect2TumorOnlyMode.err			### File to store error messages
#SBATCH --qos=normal				    	### Quality of service queue
#SBATCH --time=24:00:00						### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1 				### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwilli24@tulane.edu

# IS BENCHMARKED: NO

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

# USAGE ########################################

#

# ABOUT ########################################

# This script uses GATK Mutect2 in tumor-only mode to call variants against a reference for the purposes of making a Panel of Normals (PON) to then use in full tumor-normal somatic mutation calling.

# INPUT:
# - reference genome FASTA
# - tumor BAM

# OUTPUT:
# - raw unfiltered somatic short variant callset VCF
# - reassembled BAM
# - Mutect2 .stats file

# ARGUMENTS ########################################

in_ref_fasta=$1
in_sample_normal_bam=$2
sample_normal_bam_basename=$(basename "$in_sample_normal_bam")
out_dir=$3

# SETUP ########################################

# Make timestamped dir for results output
dir_results="${out_dir}/results--${sample_normal_bam_basename}--$(date +\%Y-\%m-\%d-%H%M%S)"
mkdir $dir_results

# Make a copy of the script run and put it into the results dir
path_script=$(echo "$0")
cp "$path_script" "$dir_results"

# MAIN CODE ########################################

## LOAD MODULES
module load gatk/4.1.8.1

## RUN ANALYSIS

# Call variants in Mutect2 tumor-only mode
gatk Mutect2 \
--reference "$in_ref_fasta" \
--input "$in_sample_normal_bam" \
--output "${dir_results}/${sample_normal_bam_basename%.bam}.normalForPON.vcf.gz"

# END ########################################
echo "End"
