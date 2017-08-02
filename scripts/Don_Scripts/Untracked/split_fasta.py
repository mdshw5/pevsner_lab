#!/usr/bin/env python3

import argparse
from io import IOBase

parser = argparse.ArgumentParser(description='Split a single fasta into seperate files')
parser.add_argument('infile', help='The input fasta file')
parser.add_argument('out_dir', help='The directory to write output files')
args = parser.parse_args()

cur_out_name = ''
cur_out = ''
cur_chrom = ''

with open(args.infile, 'r') as f:
    for line in f:
        if line[0] == '>':
            cur_chrom = line.split()[0][1:]
            cur_out_name = args.out_dir + '/' + cur_chrom + '.fa'
            if isinstance(cur_out, IOBase):
                cur_out.close()
            cur_out = open(cur_out_name, 'w')
        print(line.rstrip(), file=cur_out)
