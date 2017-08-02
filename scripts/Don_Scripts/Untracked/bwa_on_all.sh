#!/usr/bin/env bash

OUTDIR=$(pwd)

for infile in "$@"
do
    RGTWO=${infile%_1.fastq.gz}_2.fastq.gz
    BASE=${infile%_1.fastq.gz}
    BASE=${BASE##*/}

#    echo "FQ1 is $infile"
#    echo "FQ2 is $RGTWO"
#    echo "Base is $BASE"
#    echo "Output will be ${OUTDIR}/$BASE"

    bwa mem -t 20 -R '@RG\tID:foo\tSM:bar' /mnt/data/reference/hs37d5.fa $infile $RGTWO | samblaster | samtools view -Sb - > $OUTDIR/${BASE}_mem.bam
    samtools view -h $OUTDIR/${BASE}_mem.bam | samblaster -a -e --maxUnmappedBases 5 --minIndelSize 5 -d $OUTDIR/${BASE}_disc.sam -s $OUTDIR/${BASE}_sr.sam -o /dev/null
    samtools sort -@ 20 -m 3G $OUTDIR/${BASE}_mem.bam $OUTDIR/${BASE}_sorted    
done

#echo "Output directory is $OUTDIR"
