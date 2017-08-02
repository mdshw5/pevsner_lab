#!/usr/bin/env bash
#SBATCH -c 1
#SBATCH -J parse_vcf
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 7-00:05
#SBATCH -w plab-node04


for file in /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/Data_Files/VCF/*.vcf*; do
  for sample in `bcftools query -l $file`; do
    bcftools view -c1 -Oz -s $sample -o ${file/.vcf*/.$sample.vcf.gz} $file
  done
done
