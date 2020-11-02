#!/bin/bash

# use "cut" to isolate tumor VAF from a VCF file and store in a separate file for import into R to make a histogram ------------------------------

# arguments ------------------------------
VCF_FILE=$1
SAMPLE_NAME=$2


# code ------------------------------
sed '/^#/d' $VCF_FILE | cut -f 11 | cut -f 6 -d ':' | sed 's/%//g' > "$SAMPLE_NAME"-tumor-VAF.txt