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

# This script builds the ANNOVAR annotation database for the mm10 genome

# Usage ########################################

# This script should be run from within the "annovar/" folder within the experiment folder.

# experiment-folder/
# └── SAMPLE01
# └── annovar
# ├── annotate_variation.pl
# ├── coding_change.pl
# ├── convert2annovar.pl
# ├── retrieve_seq_from_fasta.pl
# ├── table_annovar.pl
# └── variants_reduction.pl
# └── build-annovar-db.sh ***

# Arguments ########################################

# Code ########################################

# 1. Build mm10 annotation database ----------

perl annotate_variation.pl \
    --downdb \
    --buildver "mm10" gene "mm10_db/"

perl annotate_variation.pl \
    --downdb \
    --buildver "mm10" seq "mm10_db/mm10_seq"

# 2. Build mm10 transcript FASTA for annotation ----------

perl retrieve_seq_from_fasta.pl \
    "mm10_db/mm10_refGene.txt" \
    --seqdir "mm10_db/mm10_seq" \
    --format refGene \
    --outfile "mm10_db/mm10_refGeneMrna.fa"

# TODO: change --format to something besides refGene (e.g. ensGene)? refGene more defined though
# OR: see
# from the <retrieve_seq_from_fasta.pl> script:
# $format =~ m/^(refGene|knownGene|ensGene|genericGene|tab|simple)$/ or pod2usage ("Error in argument: the --format argument can be only 'refGene', 'knownGene', 'ensGene', 'genericGene', 'tab' or 'simple'");

# End job ########################################
echo End Job
