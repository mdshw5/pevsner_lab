#!/usr/bin/env bash

# $1 = foo_qsorted.bam
# $2 = read group id
# $3 = sample id

name=${1%_qsorted.bam}
bedtools bamtofastq -i $1 -fq /dev/stdout /fq2 /dev/stdout | bwa mem -t 20 -R "@RG\tID:$2\tSM:$3" /mnt/data/reference/hs37d5.fa - | samblaster | samtools view -Sb - > ${name}_mem.bam
samtools view -h ${name}_mem.bam | samblaster -a -e --maxUnmappedBases 5 --minIndelSize 5 -d ${name}_disc.sam -s ${name}_sr.sam -o /dev/null
samtools view ${name}_mem.bam | tail -n+100000 | pairend_distro.py -r 100 -X 4 -N 10000 -o ${name}.pe.histo &>${name}.distro
samtools sort -@ 20 -m 1G ${name}_mem.bam ${name}_sorted_mem
samtools index ${name}_sorted_mem.bam
samtools view -Sb ${name}_disc.sam | samtools sort -@ 20 -m 1G - ${name}_disc_sorted
samtools view -Sb ${name}_sr.sam | samtools sort -@ 20 -m 1G - ${name}_sr_sorted