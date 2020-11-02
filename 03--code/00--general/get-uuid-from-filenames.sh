#!/usr/bin/env zsh

# get UUID rootnames for VCF files downloaded from GDC
cd "/Volumes/Pursell-Lab-HD/Pursell-lab/projects/project_02--DNA-seq-analysis/04--analysis/00--obtain-and-format-data/02--POLE-TCGA-cases-form-cBioPortal--Vivian/DNA-seq-Varscan2-mutations"
echo "varscan2_vcf_uuid" > ../Varscan2-vcf-uuids.tsv
for i in *.vcf ; do
	# basename=$(basename "$i")
	rootname="${i%%.*}"
	echo "$rootname" >> ../Varscan2-vcf-uuids.tsv
done

# get UUID rootname for coutns files
cd "/Volumes/Pursell-Lab-HD/Pursell-lab/projects/project_02--DNA-seq-analysis/04--analysis/00--obtain-and-format-data/02--POLE-TCGA-cases-form-cBioPortal--Vivian/RNA-seq-HTseq-counts"
echo "htseq_counts_uuid" > ../HTseq-counts-uuids.tsv
for i in *.counts ; do
	rootname="${i%%.*}"
	echo "$rootname" >> ../HTseq-counts-uuids.tsv
done