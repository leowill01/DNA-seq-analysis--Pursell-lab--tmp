#!/usr/bin/env zsh

# Description ########################################
# This script takes a folder of VCF files (e.g. multiple output VCFs from VarScan2; from mouse data) as input and uses ANNOVAR to annotate variants according to mm10 genome annotations

# Arguments ########################################

TABLE_ANNOVAR_PL_SCRIPT="table_annovar.pl" # in the ANNOVAR program folder as was downloaded
INPUT_DIR_OF_VCFS=$1
SPECIES=$2
OUTPUT_DIR=$3

if [[ "$2" == "mouse" ]]
then
	ANNOVAR_DB_DIR="mm10_db" # in the ANNOVAR program folder. made with previous 'build-annovar-db.sh' script
	ANNOVAR_REF_GENOME_SYMBOL="mm10" # e.g. 'mm10', 'hg38', etc.
elif [[ "$2" == "human" ]]
then
	ANNOVAR_DB_DIR="hg38_db" # in the ANNOVAR program folder. made with previous 'build-annovar-db.sh' script
	ANNOVAR_REF_GENOME_SYMBOL="hg38" # e.g. 'mm10', 'hg38', etc.
else
	echo "ERROR: Must enter 'mouse' or 'human' for species as second argument."
fi


# Code ########################################

# 1. Annotate mouse variants ----------

# Make results dir

	# Make timestamped run results dir inside output sample run dir
	RESULTS_DIR="${OUTPUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
	mkdir $RESULTS_DIR

# For every VCF file in `$INPUT_DIR_OF_VCFS`
for i in "${INPUT_DIR_OF_VCFS}/"**/*.vcf; do

    # extract the basename
    BASENAME=$(basename $i '.vcf')

	# Make results dir for each VCF input
	mkdir "${RESULTS_DIR}/${BASENAME}"

    # Annotate variants with ANNOVAR
    perl "$TABLE_ANNOVAR_PL_SCRIPT" \
        "$i" \
        "$ANNOVAR_DB_DIR" \
        --vcfinput \
        --outfile "${RESULTS_DIR}/${BASENAME}/${BASENAME}" \
        --buildver "$ANNOVAR_REF_GENOME_SYMBOL" \
        --protocol refGene \
        --operation g

	# Rename all .log files to .log.txt files for easier viewing
	for i in "${RESULTS_DIR}/${BASENAME}/"*".log"; do
		mv "$i" "$i".txt
	done

done

# End job ########################################

echo End Job
