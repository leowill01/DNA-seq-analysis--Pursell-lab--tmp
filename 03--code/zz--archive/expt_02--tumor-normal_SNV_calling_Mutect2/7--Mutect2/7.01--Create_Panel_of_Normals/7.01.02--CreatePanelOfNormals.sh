#!/bin/bash
#SBATCH --job-name createPoN.sh 			### Job name
#SBATCH -o Output_Log-createPoN.sh.log	### File to store output
#SBATCH -e Error_Log-createPoN.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=1					### Number of cores per node requested. Max 20/node
#SBATCH --mem=64G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Combine normal variant calls from reference into a Panel of Normals as a common germline reference for use with Mutect2


################################################################################
# ARGUMENTS
################################################################################
NORMAL_VCF_1=$1
NORMAL_VCF_2=$2
NORMAL_VCF_3=$3
NORMAL_VCF_4=$4
NORMAL_VCF_5=$5
NORMAL_VCF_6=$6
OUT_PON_PREFIX=$7


################################################################################
# EXECUTION CODE
################################################################################
module load gatk

gatk CreateSomaticPanelOfNormals \
--vcfs="$NORMAL_VCF_1" \
--vcfs="$NORMAL_VCF_2" \
--vcfs="$NORMAL_VCF_3" \
--vcfs="$NORMAL_VCF_4" \
--vcfs="$NORMAL_VCF_5" \
--vcfs="$NORMAL_VCF_6" \
--output="$OUT_PON_PREFIX".vcf \
-TMP_DIR=$TMPDIR

# This will create an 8-column *sites-only* VCF *lacking annotations*


################################################################################
# END JOB
################################################################################

echo End Job
