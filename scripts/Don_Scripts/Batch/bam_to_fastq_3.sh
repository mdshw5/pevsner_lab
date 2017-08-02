#!/usr/bin/env bash

# $1 = input
# $2 = out_base

I=$1
O=$2

#mkfifo ${O}_PIPE1
#mkfifo ${O}_PIPE2

samtools sort -n -O bam -T $O $I | bedtools bamtofastq -i /dev/stdin -fq >(gzip -c - > ${O}_1.fastq.gz)

#gzip < ${O}_PIPE1 -c > ${O}_1.fastq.gz &
#gzip < ${O}_PIPE2 -c > ${O}_2.fastq.gz

#rm ${O}_PIPE1
#rm ${O}_PIPE2
