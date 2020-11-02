#!/usr/bin/env zsh

# FILENAME:		[FILENAME]
# AUTHOR:		[AUTHOR]
# DATE:			[DATE]
# USAGE:		[USAGE]
# DESCRIPTION:	make different coverage files for testing BEDtools calculating coverage
# VERSION:		[VERSION]
# NOTES:		[NOTES]
# ZSH VERSION:	[ZSH VERSION]
# DEV PLATFORM:	[DEV PLATFORM]

# ==============================================================================

# ARGUMENTS ########################################

# BAM=$1
# BED=$2
# BED_NO_HEADER=$3
# OUT_DIR=$4

# FOR TESTING #
IN_BAM=$1
OUT_DIR=$2
BED="/Volumes/Pursell-Lab-HD/Pursell-lab/projects/project_02--DNA-seq-analysis/02--data/02--reference-data/sequencing/target-capture-regions/Agilent-SureSelect-All-Human-Exon-V6-r2/from-Agilent-website/S07604514_hs_hg38/S07604514_Regions.bed"
BED_NO_HEADER="/Volumes/Pursell-Lab-HD/Pursell-lab/projects/project_02--DNA-seq-analysis/02--data/02--reference-data/sequencing/target-capture-regions/Agilent-SureSelect-All-Human-Exon-V6-r2/from-Agilent-website/S07604514_hs_hg38/S07604514_Regions-no-header.bed.tsv"

# SETUP ########################################

# Make timestamped dir for results output
RESULTS_DIR="${OUT_DIR}/results--$(date +\%Y-\%m-\%d-%H%M%S)"
OUT_DIR_SAMT="${RESULTS_DIR}/out-samtools"
OUT_DIR_BEDT="${RESULTS_DIR}/out-bedtools"
OUT_DIR_GATK="${RESULTS_DIR}/out-gatk"
mkdir $RESULTS_DIR $OUT_DIR_SAMT $OUT_DIR_BEDT $OUT_DIR_GATK

# Make a copy of the script run and put it into the results dir
SCRIPT_PATH=$(echo "$0")
cp "$SCRIPT_PATH" "$RESULTS_DIR"

# MAIN CODE ########################################

# WORKING COMMANDS ========================================

	# SAMtools ----------------------------------------
	# Relevant commands from SAMtools manual page: http://www.htslib.org/doc/samtools.html

		# flagstat
		samtools flagstat -O tsv $IN_BAM > "${OUT_DIR_SAMT}/flagstat.tsv"

		# idxstats
		samtools idxstats $IN_BAM > "${OUT_DIR_SAMT}/idxstats.txt"

		# stats
		samtools stats "$IN_BAM" > "${OUT_DIR_SAMT}/stats.txt"
		ST_STATS="${OUT_DIR_SAMT}/stats.txt"

		samtools stats -t "$BED_NO_HEADER" "$IN_BAM" > "${OUT_DIR_SAMT}/stats-bed-regions.txt"
		ST_STATS_BED="${OUT_DIR_SAMT}/stats-bed-regions.txt"
		# TODO: test options -c, -f, -F, -t [!!: BED exome regions]

		# plot-bamstats
		plot-bamstats -p "${OUT_DIR_SAMT}/plot-bamstats/" $ST_STATS
		plot-bamstats -p "${OUT_DIR_SAMT}/plot-bamstats-bed-regions/" $ST_STATS_BED

		# bedcov: report coverage over regions in a supplied BED file
		samtools bedcov "$BED" "$IN_BAM" > "${OUT_DIR_SAMT}/bedcov.tsv"
		# only include reads with MAPQ>0
		samtools bedcov -Q 0 "$BED" "$IN_BAM" > "${OUT_DIR_SAMT}/bedcov-Qgt0.tsv"

		# depth: computes read depth at each position or region
		# relevant options: -b [BED] -f [BAMLIST.txt] -H -o [OUT.txt] -q [INT baseQ] -Q [INT MAPQ] -r [REGION CHR:FROM-TO] -g [FLAGS] -G [FLAGS]
		samtools depth -H $IN_BAM > "${OUT_DIR_SAMT}/depth.tsv" # very large filesize
		samtools depth -H -b $BED $IN_BAM > "${OUT_DIR_SAMT}/depth-b.tsv" # large filesize, but smaller bc restricted to regions in the BED file
			# use AWK to average the column of per-base reads
			echo "Average depth of coverage:" > "${OUT_DIR_SAMT}avg-depth.txt"
			gegrep '^[^#]' "${OUT_DIR_SAMT}/depth-b.tsv" | \
			awk '{total+=$3} END {print total/NR}' >> "${OUT_DIR_SAMT}avg-depth.txt"
		samtools depth -H -b $BED -Q 0 $IN_BAM > "${OUT_DIR_SAMT}/depth-b_Qgt0.tsv" # may be the same bc samtools depth automatically ignores reads with UNMAP, SECONDARY, QCFAIL, or DUP flags
		samtools depth -H -Q 0 -r "chr12:132621762-132689524" $IN_BAM > "${OUT_DIR_SAMT}/depth-rPOLE_Qgt0.tsv" # POLE gene
		samtools depth -H -Q 0 -r "chr12:132673260-132673261" $IN_BAM > "${OUT_DIR_SAMT}/depth-rPOLES459_Qgt0.tsv" # POLE:S459


	# BEDtools ----------------------------------------



	# GATK ----------------------------------------



# * NOT WORKING COMMANDS

	# # How-tos from: Kevin Blighe at https://www.biostars.org/p/279140/ 
	# # select only relevant BAM sections: https://github.com/arq5x/bedtools2/issues/492

	# # 1. get mean DOC for each interval in BED file
	# 	# TODO: can i then take a mean of those means? -> see AWK command below
	# 	# bedtools coverage -a "$bed" -b "$bam" -mean > meanCoverageBED.bedgraph &
	# 	bedtools coverage -a "$bed" -b "$bam" -mean > meanCoverageBED.txt

	# # 2. get the overall mean, pipe into AWK and get the average of the average
	# 	bedtools coverage -a "$bed" -b "$bam" -mean | awk '{total+=$4} END {print total/NR}' > overallMeanCovBED.txt

	# # 3. Get per-base read depth for each region in the BED file
	# 	# poster uses .bedgraph as filename
	# 	bedtools coverage -a "$bed" -b "$bam" -d > perBaseDepthBED.txt

# Conclusion

echo "Done"