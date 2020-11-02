#!/bin/bash

# About ########################################

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

REF_GENOME_SYMBOL="$1"

# Code ########################################

# 1. Build ref genome annotation database ----------

perl annotate_variation.pl \
    --downdb \
    --buildver "$REF_GENOME_SYMBOL" gene "${REF_GENOME_SYMBOL}_db/"

perl annotate_variation.pl \
    --downdb \
    --buildver "$REF_GENOME_SYMBOL" seq "${REF_GENOME_SYMBOL}_db/${REF_GENOME_SYMBOL}_seq"

# 2. Build ref genome transcript FASTA for annotation ----------

perl retrieve_seq_from_fasta.pl \
    "${REF_GENOME_SYMBOL}_db/${REF_GENOME_SYMBOL}_refGene.txt" \
    --seqdir "${REF_GENOME_SYMBOL}_db/${REF_GENOME_SYMBOL}_seq" \
    --format refGene \
    --outfile "${REF_GENOME_SYMBOL}_db/${REF_GENOME_SYMBOL}_refGeneMrna.fa"

# TODO: change --format to something besides refGene (e.g. ensGene)? refGene more defined though
# OR: see
# from the <retrieve_seq_from_fasta.pl> script:
# $format =~ m/^(refGene|knownGene|ensGene|genericGene|tab|simple)$/ or pod2usage ("Error in argument: the --format argument can be only 'refGene', 'knownGene', 'ensGene', 'genericGene', 'tab' or 'simple'");

# End job ########################################
echo End Job
