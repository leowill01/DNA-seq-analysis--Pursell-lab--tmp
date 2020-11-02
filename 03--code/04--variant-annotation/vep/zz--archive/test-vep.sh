# test VEP

cd ~/Desktop/tmp--work

IN_VCF="/Users/leo/Desktop/tmp--work/9238-Int.snp.Somatic.hc.filter.vcf"
OUT_TSV="${IN_VCF%.vcf}.vep.tsv"
OUT_VCF="${IN_VCF%.vcf}.vep.vcf"

printf "Arguments:\n\
IN_VCF: ${IN_VCF}\n\
OUT_TSV: ${OUT_TSV}\n\
OUT_VCF: ${OUT_VCF}"

# vep --cache --species mouse --symbol --show_ref_allele --numbers --tab -i "$IN_VCF" -o "$OUT_FILE"

# for TSV output
vep \
--verbose \
--offline \
--force_overwrite \
--cache \
--species mus_musculus \
--symbol \
--show_ref_allele \
--numbers \
--tab \
--fields "Uploaded_variation,Location,Allele,REF_ALLELE,SYMBOL,Gene,Feature,Feature_type,Consequence,IMPACT,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,Existing_variationSYMBOL_SOURCE,EXON,INTRON,DISTANCE,STRAND,FLAGS,HGNC_ID" \
-i "9238-Int.snp.Somatic.hc.filter.vcf" \
-o "9238-Int.snp.Somatic.hc.filter.vep.tsv"

# for VCF output
vep \
--verbose \
--offline \
--force_overwrite \
--cache \
--species mus_musculus \
--symbol \
--numbers \
--vcf \
-i "9238-Int.snp.Somatic.hc.filter.vcf" \
-o "9238-Int.snp.Somatic.hc.filter.vep.vcf"

# --fields ""
