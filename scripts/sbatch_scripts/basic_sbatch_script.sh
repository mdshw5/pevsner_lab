#!/bin/bash
#
#SBATCH --job-name common_experiment
#SBATCH --output /mnt/data/jeremy/10X_genomics/10X_Data/common_experiment_fastqc_output/output.txt
#SBATCH --error /mnt/data/jeremy/10X_genomics/10X_Data/common_experiment_fastqc_output/error.txt
#SBATCH --ntasks 5
#SBATCH --time 06:59:00
#SBATCH --mem-per-cpu 1000
#

#module add bwa.*
#bwa mem -t 10 /mnt/userdata/reference/hs37d5.fa /mnt/userdata/jeremy/NBU_Data/Subsets/kk_002da_subset_1.fasq /mnt/userdata/jeremy/NBU_Data/Subsets/kk_002da_subset_2.fastq \
# > /mnt/userdata/jeremy/NBU_Data/Subsets/kk002da_subset_merged_sbatch_test.sam
 

fastqc -o /mnt/data/jeremy/10X_genomics/10X_Data/common_experiment_fastqc_output /mnt/data/10X/CommonExperiment/1611UNHX-0003/PevsnerChromiumCommon1/PevsnerChromiumCommon1_R2.fastq.gz /mnt/data/10X/CommonExperiment/1611UNHX-0003/PevsnerChromiumCommon1/PevsnerChromiumCommon1_R1.fastq.gz
