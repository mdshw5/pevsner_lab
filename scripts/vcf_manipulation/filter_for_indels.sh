#!/usr/bin/env bash
#SBATCH -c 1
#SBATCH -J parse_vcf
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 7-00:05
#SBATCH -w plab-node04


for file in /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/Data_Files/VCF/*.vcf.gz;
do
    echo $file
    basename=${file##*combined_raw.sorted.}
    echo $basename
    vcfname=${basename%%.gz}
    outname=${basename%%.vcf.gz}
    echo $outname
    vcftools --gzvcf $file --out /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/Data_Files/VCF/$outname.indels --keep-only-indels --recode
    bgzip /mnt/data/jeremy/projects/sws/collab_samples/gatk_output/Data_Files/VCF/$outname.indels.recode.vcf
done
