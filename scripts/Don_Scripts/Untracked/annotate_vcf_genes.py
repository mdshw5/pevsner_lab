#!/usr/bin/env python3

import argparse
import sys

parser = argparse.ArgumentParser(description='Annotate genes marked in a VCF based on if they can be found in a list of known candidates.')
parser.add_argument('gene_list', help='A list of genes to cross-reference. One gene per line.')
parser.add_argument('tag', help='A tag to put into the INFO field of the VCF.')
parser.add_argument('--input', help='An input file (default is stdin)')
args = parser.parse_args()

genes = {}

with open(args.gene_list, 'r') as f:
    for line in f:
        line = line.rstrip()
        genes[line] = True

if args.input:
    infile = open(args.input, 'r')
else:
    infile = sys.stdin

for line in infile:
    line = line.rstrip().split()
    if 'MISSENSE' in line[7]:
        gene_idx = line[7].index('MISSENSE')
        idx_end = line[7].index(';', gene_idx)
        gene = line[7][gene_idx + 14 : idx_end]
        if gene in genes:
            line[7] = args.tag + ';' + line[7]
    elif 'INTRONIC' in line[7]:
        gene_idx = line[7].index('INTRONIC')
        idx_end = line[7].index(';', gene_idx)
        gene = line[7][gene_idx + 9 : idx_end]
        if gene in genes:
            line[7] = args.tag + ';' + line[7]
    elif 'NONSENSE' in line[7]:
        gene_idx = line[7].index('NONSENSE')
        idx_end = line[7].index(';', gene_idx)
        gene = line[7][gene_idx + 9 : idx_end]
        if gene in genes:
            line[7] = args.tag + ';' + line[7]
    print('\t'.join(line))
