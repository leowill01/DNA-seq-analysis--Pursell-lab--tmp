#!/bin/bash

echo Start Job

################################################################################
# JOB DESCRIPTION
################################################################################
# Uses ANNOVAR to build a nonhuman transcript FASTA file for variant annotation.


################################################################################
# ARGUMENTS
################################################################################


################################################################################
# EXECUTION CODE
################################################################################
perl retrieve_seq_from_fasta.pl \
"mm10db"/"mm10"_refGene.txt \
--seqdir "mm10db"/"mm10"_seq \
--format refGene \
--outfile "mm10db"/"mm10"_refGeneMrna.fa

# TODO: change --format to ensGene instead of refGene? refGene more defined though
    # OR: see 
# from the <retrieve_seq_from_fasta.pl> script:
    # $format =~ m/^(refGene|knownGene|ensGene|genericGene|tab|simple)$/ or pod2usage ("Error in argument: the --format argument can be only 'refGene', 'knownGene', 'ensGene', 'genericGene', 'tab' or 'simple'");


################################################################################
# END JOB
################################################################################
echo End Job
