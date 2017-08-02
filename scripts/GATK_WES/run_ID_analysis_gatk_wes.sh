#!/usr/bin/env bash

module load anaconda3/2.3.0 anaconda bwa samtools sambamba samblaster GATK JDK tabix_bgzip bcftools svtyper lumpy/0.2.13 freebayes

# Document
dt=$(sed 's/-//g'<<< $(date "+%F [%T]"))
echo "${dt} : Begin GATK WES analysis [${0}]" >> CHANGELOG

# Run GATK

snakemake -w 60 --configfile Configs/ID_gatk.json -p -j 30 --nocolor --verbose -s pipeline/gatk_wes_call_recal_vars.smk --cluster "sbatch -N 1 -n 1 -c {threads} --error Logs/gatk_{rule}_%j.log --output Logs/gatk/gatk_{rule}_%j.log"


#Removed sbatch "nodelist=plabnode04" from snakemake command
