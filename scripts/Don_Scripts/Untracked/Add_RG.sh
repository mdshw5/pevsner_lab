#!/usr/bin/env bash

for infile in "$@"
do
    java -jar -Xmx4g $PICARD/AddOrReplaceReadGroups.jar INPUT\= $infile OUTPUT\= ${infile%_sorted_bwa.bam}.rg.bam RGLB\=Control_Exome RGPL\=illumina RGPU\=blue RGSM\=${infile%_sorted_bwa.bam}
done