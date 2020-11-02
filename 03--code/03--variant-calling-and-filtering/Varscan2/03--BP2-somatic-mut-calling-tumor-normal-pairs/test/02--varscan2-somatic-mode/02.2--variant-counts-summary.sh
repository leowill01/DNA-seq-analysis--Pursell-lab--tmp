#!/bin/bash

echo 'Total SNP Counts' > summary-total-counts.txt

for i in */*.snp.vcf
do 
    echo $i | cut -f 1 -d '/' >> summary-total-counts.txt
    sed '/^#/d' $i | wc -l >> summary-total-counts.txt
done

echo 'Total Indel Counts' >> summary-total-counts.txt

for i in */*.indel.vcf
do 
    echo $i | cut -f 1 -d '/' >> summary-total-counts.txt
    sed '/^#/d' $i | wc -l >> summary-total-counts.txt
done
