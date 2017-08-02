#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser(description='Quickly find if variants overlap regions of interest')
parser.add_argument('GFF_file', help='A file in GFF3 format with exons and genes')
parser.add_argument('chromosome', help='Variant chromosome')
parser.add_argument('position', help='Variant position')
args = parser.parse_args()

chroms = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'MT':25, 'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24, 'chrM':25}

with open(args.GFF_file, 'r') as f:
    for line in f:
        if line[0] == '#':
            continue
        line = line.rstrip().split()
        if not line[0] in chroms:
            sys.exit('Chromosome not found')
        if chroms[line[0]] == chroms[args.chromosome]:
            if int(args.position) < int(line[3]):
                sys.exit(0)
            else:
                if int(args.position) < int(line[4]):
                    print('\t'.join(line))
                else:
                    pass
        else:
            continue
