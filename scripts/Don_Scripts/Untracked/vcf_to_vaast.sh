#!/usr/bin/env bash

#echo "$0 $1 $2 $3"

vaast_converter -n 0 --build hg19 --path $2 $1

count=0
files=''

for infile in `ls $2*.gvf`
do
    base=${infile%.gvf}
    /mnt/data/don/Scripts/gffsort.pl $infile > ${base}.sorted.gvf
    VAT --build hg19 -f /mnt/data/reference/VAAST/refGene_hg19_with_introns.header.sorted.gff3 -a /mnt/data/reference/hg19/ucsc.hg19.fasta ${base}.sorted.gvf > ${base}.vat.gvf
    files="${base}.vat.gvf $files"
    count=$(($count + 1))
done

count=$(($count - 1))

VST -o "U(0..${count})" -b hg19 $files > $3

#VAAST -gp 1e6 -iht n -m lrt -o vaast_out -p 12 --use_aas_info y --splice_site -r 0.05 -k /mnt/data/reference/VAAST/refGene_hg19_with_introns.header.sorted.gff3 /mnt/data/reference/VAAST/1KG_refGene_Dec2011_CGDiv_NHLBI_NoCall.cdr all_files.cdr

