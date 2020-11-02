#!/bin/bash
#SBATCH --job-name analyzeCovariates.sh 			### Job name
#SBATCH -o Output_Log-analyzeCovariates.sh.log	### File to store output
#SBATCH -e Error_Log-analyzeCovariates.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE benchmarked


echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Uses GATK AnalyzeCovariates to visualize the quality of a recalibration run done by GATK BaseRecalibrator. It uses at minimum a recalibration table from the original alignment (i.e. before the application of recalibrationTable1.table to the original BAM file with ApplyBQSR) and a 2nd recalibration table after applying recalibrationTable1.table to the original BAM file.

################################################################################
# ARGUMENTS
################################################################################
RECALIBRATION_TABLE_1=$1
RECALIBRATION_TABLE_2=$2
OUT_CSV_AND_PDF_HEADER=$3

################################################################################
# EXECUTION CODE
################################################################################
module load gatk
module load R

gatk AnalyzeCovariates \
-before "$RECALIBRATION_TABLE_1" \
-after "$RECALIBRATION_TABLE_2" \
-csv "$OUT_CSV_AND_PDF_HEADER".BQSR_R1-R2.csv \
-plots "$OUT_CSV_AND_PDF_HEADER".BQSR_R1-R2.pdf \
-TMP_DIR=$TMPDIR

################################################################################
# END JOB
################################################################################
echo End Job
