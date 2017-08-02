#!/usr/bin/env bash

# $1 = foo_qsorted.bam
# $2 = read group id
# $3 = sample id
# $4 = directory
# $5 = cores

name=${1%.bam}
name=${name##*/}
samtools sort -n -@ $5 -m 1G $1 $4/${name}_qsorted
bedtools bamtofastq -i $4/${name}_qsorted.bam -fq /dev/stdout -fq2 /dev/stdout | bwa mem -p -t $5 -R "@RG\tID:$2\tSM:$3" /mnt/data/reference/hs37d5.fa - | samblaster | samtools view -Sb - > $4/${name}_mem.bam
rm $4/${name}_qsorted.bam
samtools view -h $4/${name}_mem.bam | samblaster -a -e --maxUnmappedBases 5 --minIndelSize 5 -d $4/${name}_disc.sam -s $4/${name}_sr.sam -o /dev/null
samtools view $4/${name}_mem.bam | tail -n+100000 | pairend_distro.py -r 100 -X 4 -N 10000 -o $4/${name}.pe.histo &>$4/${name}.distro
samtools sort -@ $5 -m 1G $4/${name}_mem.bam $4/${name}_sorted_mem
rm $4/${name}_mem.bam
samtools index $4/${name}_sorted_mem.bam
samtools view -Sb $4/${name}_disc.sam | samtools sort -@ $5 -m 1G - $4/${name}_disc_sorted
samtools view -Sb $4/${name}_sr.sam | samtools sort -@ $5 -m 1G - $4/${name}_sr_sorted