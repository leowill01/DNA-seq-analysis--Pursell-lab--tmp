#!/bin/bash
#SBATCH --job-name=trint	                ### Job name
#SBATCH --output=Output_Log-trint.log	### File to store output
#SBATCH --error=Error_Log-trint.err	    ### File to store error messages
#SBATCH --qos=normal				    ### Quality of service queue
#SBATCH --time=24:00:00					### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1						### Nodes requested for job
#SBATCH --ntasks-per-node=1 			### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1				### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G						### Memory requested for job. Max 64GB/node
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lwilli24@tulane.edu

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

# About ##############################
# This script take a vcf input and adds a column to each variant with the trinucleotide context

# NOTE: Requires Bash 4.x+ to run due to usage of `shopt -s globstar`
# NOTE: Requires Python 3.x+ to run Nate's trinucleotide-context.py script

# Usage ##############################
# Where: working dir = ".../expt_05--trinucleotide-context/sample_01/"
# How: bash script.sh [IN_VCF_DIR/] &> /dev/null &

# Arguments ##############################
CHR_FASTA_DIR="/lustre/project/zpursell/leo/projects/project-02--DNA-seq-analysis-of-mouse-POLE-exo-mutant-tumors/02--raw-data/02--reference-data/genome/mm10--UCSC/mm10_chr_fa_UCSC"
IN_VCF_DIR=$1

# Code ##############################
# Load Anaconda module on Cypress to access Python 3.x+
# NOTE: Tulane HPC says to use `export CONDA_ENVS_PATH=/lustre/project/zpursell/leo/conda-envs` to store envs in project folder (I guess for easier sharing among labmates?), but for now I'm just leaving them in the default `/home/lwilli24/.conda/envs/[env-name]` dir.
module load anaconda3/5.1.0

# Switch to conda environment that has bioinformatics dependencies (e.g. Biopython) installed (as per Tulane HPC wiki "Install Packages to Anaconda Python"). Whatever environment you're switching to must already exist.
source activate env-bioinformatics

# Enable globstar in Bash (v4.x+)
shopt -s globstar

# Make dir to store results
mkdir results

# Run trinucleotide python script for every VCF file in a directory containing output VCF files (e.g. from a variant calling or VCF annotation experiment)
for i in ${IN_VCF_DIR}/**/*snp*.vcf; do
	BASE=$(basename "$i")
	python ../trinucleotide-context.py $i $CHR_FASTA_DIR > "results/${BASE}.trint.txt"
done