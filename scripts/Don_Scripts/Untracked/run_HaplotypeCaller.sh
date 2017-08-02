#!/usr/bin/env bash

INPUT=''

for infile in "$@"
do
    samtools index $infile
    INPUT="${INPUT}-I $infile "
done

java -Xmx100g -jar $GATK -T HaplotypeCaller -R /mnt/data/reference/hs37d5.fa -nct 25 --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 20 -D /mnt/data/reference/GATK_b37/dbsnp_137.b37.vcf -o haplotypeCaller_out.vcf -A BaseQualityRankSumTest -A RMSMappingQuality -A TandemRepeatAnnotator -A QualByDepth -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A Coverage -A MappingQualityZero -A SpanningDeletions $INPUT

/mnt/data/don/Scripts/filter_vcf.py -BQRS 4 -depth 660 -FSBias 4 -mapq 10 -MQRS 4 -ReadPosRS 4 -QD 10 haplotypeCaller_out.vcf variants_annotated.vcf

grep '^#\|PASS' variants_annotated.vcf > variants_filtered.vcf