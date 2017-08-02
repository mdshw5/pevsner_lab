#!/usr/bin/env python3

import sys
import argparse

parser = argparse.ArgumentParser(description='Remove variants that are not exonic.')
parser.add_argument('gff3_file', help='A refence file containing exons')
args = parser.parse_args()

ref_file = open(args.gff3_file, 'r')

chroms = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'MT':25, 'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24, 'chrM':25}

for line in sys.stdin:
    line = line.rstrip().split()
    if not line[0] in chroms:
        continue
    ref_file.seek(0)
    for ref_line in ref_file:
        if ref_line[0] == '#':
            continue
        ref_line = ref_line.rstrip().split()
        if ref_line[2] != 'exon':
            continue
        if not ref_line[0] in chroms:
            break
        if chroms[line[0]] == chroms[ref_line[0]]:
            if int(line[1]) < int(ref_line[3]):
                break
            else:
                if int(line[1]) < int(ref_line[4]):
                    print('\t'.join(line))
                    break
        else:
            continue
