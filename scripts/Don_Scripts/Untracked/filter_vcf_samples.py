#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(description='Output a vcf containing only samples found in another file')
parser.add_argument('sample_file', help='A file with samples to filter on, one sample per line')
parser.add_argument('vcf_file', help='The vcf file to filter.')
parser.add_argument('output', help='An output vcf file')
args = parser.parse_args()

sample_dict = {}

with open(args.sample_file, 'r') as f:
    for line in f:
        line = line.rstrip()
        sample_dict[line] = 0

outfile = open(args.output, 'w')

output_fields = [0, 1, 2, 3, 4, 5, 6, 7, 8]

with open(args.vcf_file, 'r') as f:
    for line in f:
        line = line.rstrip()
        if line[0:2] == '##':
            print(line, file=outfile)
            continue
        line = line.split()
        if line[0] == '#CHROM':
            for i in range(len(line)):
                if line[i] in sample_dict:
                    output_fields.append(i)
        out_array = [ line[x] for x in output_fields ]
        print('\t'.join(out_array), file=outfile)
