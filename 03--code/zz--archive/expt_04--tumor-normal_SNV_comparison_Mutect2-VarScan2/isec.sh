#!/bin/bash

# ARGUMENTS
ISEC_FOLDER_LABEL=$1


# CODE
for i in *.vcf ; do
	bgzip $i
done

for v in *.vcf.gz ; do
	tabix $v
done

bcftools isec \
*mutect*.gz \
*varscan*.gz \
-p "$ISEC_FOLDER_LABEL"
