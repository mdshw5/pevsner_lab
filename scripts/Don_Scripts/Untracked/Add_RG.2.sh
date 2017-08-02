#!/usr/bin/env bash

for infile in "$@"
do
    java -jar -Xmx4g $PICARD/AddOrReplaceReadGroups.jar INPUT\= $infile OUTPUT\= ${infile%.bam}_rg.bam RGLB\=Public_genomes RGPL\=illumina RGPU\=blue RGSM\=${infile%.bam}
done