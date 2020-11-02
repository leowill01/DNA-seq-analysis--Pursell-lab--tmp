#!/usr/bin/env zsh

# title: 01--annotate-variants-vep.sh
# author: Leo Williams
# date: 2020-06-14

# USAGE ########################################

	# zsh 01--annotate-variants-vep.sh $in_dirh_vcfs $vep_species $ref_fasta

# ABOUT ########################################

	# This script takes a dir of VCF files and annotates them with the Variant Effect Predictor (VEP) command line tool and outputs an annotated VCF as well as a human-readable TSV. For info on options, refer to: http://useast.ensembl.org/info/docs/tools/vep/script/vep_options.html#opt_af.

	# Options to explore:
	# --everything \
	# --hgvs \
	# --hgvsg \
	# --transcript_version \
	# --protein \
	# --ccds \
	# --uniprot \
	# --domains \
	# --no_check_alleles \
	# --af \
	# --max_af \
	# --pubmed \

# ARGUMENTS ########################################

	in_dirh_vcfs=$1 # RECURSIVE, annotates VCF files and outputs to same dir as VCF with added `*.vep.vcf` extension
	species_common_name=$2 # either 'mouse' or 'human'
	out_dir=$3

	# select correct reference genome for entered species
	if [[ "$species_common_name" = "mouse" ]] ; then
		vep_species="mus_musculus"
		ref_fasta="/Volumes/Pursell-Lab-HD/Pursell-lab/02--projects/project_02--DNA-seq-analysis/02--data/02--reference-data/genome/mm10--UCSC/mm10_UCSC.fa"
	elif [[ "$species_common_name" = "human" ]] ; then
		vep_species="homo_sapiens"
		ref_fasta="/Volumes/Pursell-Lab-HD/Pursell-lab/02--projects/project_02--DNA-seq-analysis/02--data/02--reference-data/genome/hg38--UCSC/rsync/hg38--UCSC/hg38.fa"
	fi

	printf "Arguments:\n\
	in_dirh_vcfs=${in_dirh_vcfs}\n\
	vep_species=${vep_species}\n\
	ref_fasta=${ref_fasta}\n"

# SETUP ########################################

	# Make timestamped dir for results output
	dir_results="${out_dir}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
	mkdir $dir_results

	# Make a copy of the script run and put it into the results dir
	path_script=$(echo "$0")
	cp "$path_script" "$dir_results"

# MAIN CODE ########################################

	# For each VCF file recursively within the directory, annotate with VEP and output an annotated VCF.vep.vcf and easily-readable TSV.vep.tsv
	for ivcf in $in_dirh_vcfs/**/*.vcf ; do
		# get basename for output file
		ivcf_basename=$(basename "$ivcf")

		# specify which file is being analyzed
		printf "\nAnalyzing file:\n\
		${ivcf}\n\n"

		# Annotate and output a human-readable TSV with all features
		vep \
		--verbose \
		--offline \
		--force_overwrite \
		--cache \
		--format vcf \
		--species "$vep_species" \
		--symbol \
		--show_ref_allele \
		--flag_pick \
		--fields "Uploaded_variation,Location,REF_ALLELE,Allele,SYMBOL,Protein_position,Consequence,IMPACT,EXON,INTRON,Gene,Feature,Feature_type,cDNA_position,CDS_position,Amino_acids,Codons,Existing_variation,DISTANCE,STRAND,FLAGS,SYMBOL_SOURCE,HGNC_ID,CCDS,ENSP,SWISSPROT,TREMBL,UNIPARC,DOMAINS,HGVSc,HGVSp,HGVS_OFFSET,HGVSg,AF,MAX_AF,MAX_AF_POPS,CLIN_SIG,SOMATIC,PHENO,PUBMED" \
		--numbers \
		--tab \
		--no_stats \
		--hgvs \
		--hgvsg \
		--transcript_version \
		--protein \
		--ccds \
		--uniprot \
		--domains \
		--no_check_alleles \
		--af \
		--max_af \
		--pubmed \
		--fasta "$ref_fasta" \
		-i "$ivcf" \
		-o "$dir_results"/"${ivcf_basename%.vcf}.vep.all.tsv"

		# Annotate and output a human-readable TSV with single annotation per variant based on rank order criteria found at http://useast.ensembl.org/info/docs/tools/vep/script/vep_other.html#summary_eg. Primary criteria is canonical transcript status.
		vep \
		--verbose \
		--offline \
		--force_overwrite \
		--cache \
		--format vcf \
		--species "$vep_species" \
		--symbol \
		--show_ref_allele \
		--pick \
		--fields "Uploaded_variation,Location,REF_ALLELE,Allele,SYMBOL,Protein_position,Consequence,IMPACT,EXON,INTRON,Gene,Feature,Feature_type,cDNA_position,CDS_position,Amino_acids,Codons,Existing_variation,DISTANCE,STRAND,FLAGS,SYMBOL_SOURCE,HGNC_ID,CCDS,ENSP,SWISSPROT,TREMBL,UNIPARC,DOMAINS,HGVSc,HGVSp,HGVS_OFFSET,HGVSg,AF,MAX_AF,MAX_AF_POPS,CLIN_SIG,SOMATIC,PHENO,PUBMED" \
		--numbers \
		--tab \
		--no_stats \
		--hgvs \
		--hgvsg \
		--transcript_version \
		--protein \
		--ccds \
		--uniprot \
		--domains \
		--no_check_alleles \
		--af \
		--max_af \
		--pubmed \
		--fasta "$ref_fasta" \
		-i "$ivcf" \
		-o "$dir_results"/"${ivcf_basename%.vcf}.vep.pick.tsv"

		# for VCF output
		vep \
		--verbose \
		--offline \
		--force_overwrite \
		--cache \
		--format vcf \
		--flag_pick \
		--show_ref_allele \
		--fields "Uploaded_variation,Location,REF_ALLELE,Allele,SYMBOL,Protein_position,Consequence,IMPACT,EXON,INTRON,Gene,Feature,Feature_type,cDNA_position,CDS_position,Amino_acids,Codons,Existing_variation,DISTANCE,STRAND,FLAGS,SYMBOL_SOURCE,HGNC_ID,CCDS,ENSP,SWISSPROT,TREMBL,UNIPARC,DOMAINS,HGVSc,HGVSp,HGVS_OFFSET,HGVSg,AF,MAX_AF,MAX_AF_POPS,CLIN_SIG,SOMATIC,PHENO,PUBMED" \
		--species "$vep_species" \
		--symbol \
		--numbers \
		--vcf \
		--hgvs \
		--hgvsg \
		--transcript_version \
		--protein \
		--ccds \
		--uniprot \
		--domains \
		--no_check_alleles \
		--af \
		--max_af \
		--pubmed \
		--fasta "$ref_fasta" \
		-i "$ivcf" \
		-o "$dir_results"/"${ivcf_basename%.vcf}.vep.vcf"

		# --fields ""

	done

printf "\nDone\n"

	# Q
	# --flag_pick for VCF?
