#!/usr/bin/env bash

#module load anaconda3 anaconda bwa samtools sambamba samblaster GATK JDK tabix_bgzip bcftools svtyper lumpy/0.2.13 freebayes

# GATK config
#python3 pipeline/make_unit_config.py --fastq data/fastq_processed/*R1.fastq.gz --config  Configs/ID_gatk.json GATK $GATK rm_rep.py

python3 make_unit_config.py --fastq /mnt/data/jeremy/projects/sws/collab_samples/data/*R1.fastq.gz --config  ID_gatk.json GATK $GATK rm_rep.py
# Lumpy config

#python3 pipeline/make_config.py --fastq data/fastq_processed/*R1.fastq.gz --config Configs/ID_lumpy.json lumpy `which svtyper` pipeline/lumpy_exclude.py

# Freebayes config

#python3 pipeline/make_config.py --fastq data/fastq_processed/*R1.fastq.gz --config Configs/ID_freebayes.json freebayes pipeline/rm_rep.py

# Document
#dt=$(sed 's/-//g'<<< $(date "+%F [%T]"))
#echo "${dt} : Make snakemake Configs [${0}]" >> CHANGELOG
