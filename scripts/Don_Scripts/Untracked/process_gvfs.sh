#!/usr/bin/env bash

base=${1%.gvf}
/mnt/data/don/Scripts/gffsort.pl $1 > ${base}.sorted.gvf
VAT --build hg19 -f /mnt/data/reference/VAAST/refGene_hg19_with_introns.header.sorted.gff3 -a /mnt/data/reference/hg19/ucsc.hg19.fasta ${base}.sorted.gvf > ${base}.vat.gvf