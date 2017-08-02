#!/usr/bin/env bash

#Requires:
#samtools/0.1.19
#bedtools/2.20.1
#bwa/0.7.9a
#samblaster/0.1.20
#lumpy/0.2.5
#JDK/1.8.0_11
#GATK/3.2-2
#python2_lib/0.1.0

####  PLEASE CHANGE THE READ LENGTH IN LINE 30 BEFORE USE  ######################################

INPUT=''

count=0
for infile in "$@"
do
    count=$(($count + 1))
    base=${infile%.bam}
    DIR="`pwd`/"
    base=${base##*/}
    samtools sort -n -@ 8 -m 1G $infile ${DIR}${base}_qsorted
    bedtools bamtofastq -i ${DIR}${base}_qsorted.bam -fq /dev/stdout -fq2 /dev/stdout | bwa mem /mnt/userdata/reference/hs37d5.fa -R "@RG\tID:foo\tSM:$base$count" -t 8 -p - | /data/don/PROGRAMS/samblaster/samblaster | samtools view -Sb - > ${DIR}${base}_mem.bam
    samtools view -h ${DIR}${base}_mem.bam | /data/don/PROGRAMS/samblaster/samblaster -a -e --maxUnmappedBases 5 --minIndelSize 5 -d ${DIR}${base}_disc.sam -s ${DIR}${base}_sr.sam -o /dev/null
    samtools sort -@ 8 -m 1G ${DIR}${base}_mem.bam ${DIR}${base}_sorted_mem
    samtools index ${DIR}${base}_sorted_mem.bam
    samtools view -Sb ${DIR}${base}_disc.sam | samtools sort -@ 8 -m 1G - ${DIR}${base}_disc_sorted
    samtools view -Sb ${DIR}${base}_sr.sam | samtools sort -@ 8 -m 1G - ${DIR}${base}_sr_sorted
    samtools view ${DIR}${base}_mem.bam | tail -n+100000 | /data/don/PROGRAMS/Lumpy/lumpy-sv/scripts/pairend_distro.pl -rl 100 -X 4 -N 10000 -o ${DIR}${base}.pe.histo &>${DIR}${base}.distro
    INPUT="${INPUT}-I ${DIR}${base}_sorted_mem.bam "
done

java -Xmx27g -jar /data/don/PROGRAMS/GenomeAnalysisTK-3.2.2/GenomeAnalysisTK.jar -T HaplotypeCaller -R /mnt/userdata/reference/hs37d5.fa -nct 8 --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 20 -D /mnt/userdata/reference/GATK_b37/dbsnp_137.b37.vcf -o variants_out.vcf -A BaseQualityRankSumTest -A RMSMappingQuality -A TandemRepeatAnnotator -A QualByDepth -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A Coverage -A MappingQualityZero -A SpanningDeletions $INPUT

#java -Xmx10g -jar $GATK -R /mnt/data/reference/hs37d5.fa -T VariantFiltration -o variants_annotated.vcf --variant variants_filtered.vcf --filterName "BASEQRS" --filterExpression "BaseQRankSum > 4.0 || BaseQRankSum < -4.0" --filterName "CLIPRS" --filterExpression "ClippingRankSum > 4.0 || ClippingRankSum < -4.0" --filterName "DEPTH" --filterExpression "DP < 121" --filterName "FSBIAS" --filterExpression "FS > 3.0 || FS < -3.0" --filterName "MAPQ" --filterExpression "MQ < 6" --filterName "MAPQRS" --filterExpression "MQRankSum < -3.0 || MQRankSum > 3.0" --filterName "QUALD" --filterExpression "QD < 10" --filterName "READPOSRS" --filterExpression "ReadPosRankSum < -3.0 || ReadPosRankSum > 3.0"

/mnt/userdata/don/Scripts/filter_vcf.py -BQRS 4 -depth 305 -FSBias 4 -mapq 10 -MQRS 4 -ReadPosRS 4 -QD 10 variants_out.vcf variants_annotated.vcf

grep '^#\|PASS' variants_annotated.vcf > variants_filtered.vcf