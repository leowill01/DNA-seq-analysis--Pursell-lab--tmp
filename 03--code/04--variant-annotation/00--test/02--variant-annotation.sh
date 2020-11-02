#!/bin/bash
#SBATCH --job-name=var-anno	                ### Job name
#SBATCH --output=Output_Log-var-anno.log	### File to store output
#SBATCH --error=Error_Log-var-anno.err	    ### File to store error messages
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

# Description ########################################
# This script takes a folder of VCF files (e.g. multiple output VCFs from VarScan2; from mouse data) as input and uses ANNOVAR to annotate variants according to mm10 genome annotations

# Usage ########################################
# This script should be run from within the sample-specific folder within the experiment folder:

# experiment-folder/
# └── SAMPLE01
    # └── annotate-variants.sh <<<
# └── annovar
    # ├── annotate_variation.pl
    # ├── coding_change.pl
    # ├── convert2annovar.pl
    # ├── retrieve_seq_from_fasta.pl
    # ├── table_annovar.pl
    # └── variants_reduction.pl

# Arguments ########################################

IN_FOLDER_OF_VCFS=$1

# Code ########################################

# 1. Annotate mouse variants ----------

# Make dir to store results
mkdir results

# For every VCF file in `$IN_FOLDER_OF_VCFS`
for i in $IN_FOLDER_OF_VCFS/*.vcf; do

    # extract the basename
    BASE=$(basename $i)

	# Extract file root
	ROOT=${BASE%.vcf}

	# Make results dir for each VCF input
	mkdir "results/${ROOT}"

    # Annotate variants with ANNOVAR
    perl ../annovar/table_annovar.pl \
        "$i" \
        "../annovar/mm10_db/" \
        --vcfinput \
        --outfile "results/${ROOT}/${ROOT}" \
        --buildver "mm10" \
        --protocol refGene \
        --operation g > "results/${ROOT}/${ROOT}.log.txt"
done

# End job ########################################

echo End Job