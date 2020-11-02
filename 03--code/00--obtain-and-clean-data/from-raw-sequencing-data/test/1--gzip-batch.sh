#!/bin/bash
#SBATCH --job-name gzip.sh 			### Job name
#SBATCH -o Output_Log-gzip.sh.log	### File to store output
#SBATCH -e Error_Log-gzip.sh.err	### File to store error messages
#SBATCH --qos normal						### Quality of service queue
#SBATCH -t 24:00:00							### Time limit in [days]-[hh]:[mm]:[ss]
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
# Uses gzip to unzip multiple .fq.gz files in the same folder into .fq files

################################################################################
# ARGUMENTS
################################################################################
# --

################################################################################
# EXECUTION CODE
################################################################################
for i in *.fq.gz; do
	gzip -d < $i > $(echo $i | cut -f 1 -d ".").fq
done

################################################################################
# END JOB
################################################################################
echo End Job
