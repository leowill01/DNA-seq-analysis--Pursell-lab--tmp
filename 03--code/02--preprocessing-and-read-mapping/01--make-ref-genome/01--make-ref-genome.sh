#!/bin/bash
#SBATCH --job-name makeRefGenome.sh 			### Job name
#SBATCH -o Output_Log-makeRefGenome.sh.log	### File to store output
#SBATCH -e Error_Log-makeRefGenome.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 12:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
#SBATCH --nodes=1							### Nodes requested for job
#SBATCH --ntasks-per-node=1					### Number of tasks to be launched per node
#SBATCH --cpus-per-task=10					### Number of cores per node requested. Max 20/node
#SBATCH --mem=5G							### Memory requested for job. Max 64GB/node
#SBATCH --mail-type ALL
#SBATCH --mail-user lwilli24@tulane.edu

# NOTE: Benchmark: 10 CPU cores, 5GB memory (50min) for BWA index only

echo Start Job
echo nodes: $SLURM_JOB_NODELIST
echo job id: $SLURM_JOB_ID
echo Number of tasks: $SLURM_NTASKS

################################################################################
# JOB DESCRIPTION
################################################################################
# Uses BWA to index reference genome FASTA file.
# Usage: sbatch <script> <genome.fa>

# Creates a dictionary for the reference genome for downstream use of other Picard & GATK4 tools.

# Make an index for the reference genome FASTA with SAMtools. This is needed for use by various GATK tools such as HaplotypeCaller


################################################################################
# ARGUMENTS
################################################################################
REF_FASTA=$1
GATK4_EXEC="/lustre/project/zpursell/leo/code-packages-external/gatk-4.1.2.0/gatk"


################################################################################
# EXECUTION CODE
################################################################################
module load bwa
module load java-openjdk/1.8.0 # load Java version required for recent Picard version
module load samtools # must use Cypress version bc macOS-built version doesnt work

bwa index \
-a bwtsw \
"$REF_FASTA"

# index ref genome FASTA with SAMtools
samtools faidx "$REF_FASTA"

# Make basename for dict file (need to remove the ".fa" when adding ".dict" because GATK tools won't recognize "ref.fa.dict", only "ref.dict")
BASE="${REF_FASTA%.*}"
echo $BASE

module load java-openjdk/1.8.0 # Java version needed for newest Picard version
# make ref genome FASTA dictionary with Picard
$GATK4_EXEC CreateSequenceDictionary \
-R "$REF_FASTA" \
-O "${BASE}.dict"


################################################################################
# END JOB
################################################################################
echo End Job
