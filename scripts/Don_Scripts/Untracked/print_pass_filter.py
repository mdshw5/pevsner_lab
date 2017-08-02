#!/usr/bin/env python3

import argparse
import shlex
import subprocess
import sys
import os.path
import time

parser = argparse.ArgumentParser(description='Takes a list of bam files as input, realigns them using bwa mem, and calls variants using LUMPY.')
parser.add_argument('--filter_file', help='A file used to define which inputs to process')
parser.add_argument('--filter_sample_column', help='In the filter file, a column which as the sample id. Sample id must = input file base name.')
parser.add_argument('--filter_text', nargs='+', help='Text which must be found in the sample line of the filter file for that sample to be output.')
parser.add_argument('inputs', nargs='+', help='Input files in .bam format.')
args = parser.parse_args()

samples = {}

if args.filter_file:
    if not args.filter_sample_column or not args.filter_text:
        sys.exit("If a filtering file is supplied --filter_sample_column and --filter_text are required arguments")
    with open(args.filter_file, 'r') as f:
        for line in f:
            print_line = [False] * len(args.filter_text)
            line = line.rstrip().split('\t')
            for column in line:
                for i in range(len(args.filter_text)):
                    if column == args.filter_text[i]:
                        print_line[i] = True
            samples[line[int(args.filter_sample_column)]] = all(print_line)

for infile in args.inputs:
    if infile.endswith('.bam'):
        base = infile[infile.rfind('/') + 1 : -4]
    else:
        print('Input file does not end with ".bam". Is the input file in bam format? Skipping: ' + infile, file=sys.stdout)
        continue
    if not base in samples:
        print(str(base) + ' not found in filtering file', file=sys.stderr)
    else:
        if samples[base]:
            print(base)
