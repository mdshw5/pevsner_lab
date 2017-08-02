#!/usr/bin/env python

from __future__ import print_function

import argparse

parser = argparse.ArgumentParser(description='Find sra files which are not unpacked')
parser.add_argument('-infile', help='A list of all of the input .sra files separated by line.')
parser.add_argument('-outfile', help='A file to write the list of unpacked files to.')
parser.add_argument('-bam_list', help='A list of output bam files, seperated by line.')
args = parser.parse_args()

bam_dict = {}

with open(args.bam_list, 'r') as f:
    for line in f:
        line = line.rstrip()
        if line[-4:] == '.bam':
            m_line = line[-13:]
            print(m_line[:-4])
            bam_dict[m_line[:-4]] = 0

print(str(bam_dict))

outfile = open(args.outfile, 'w')
with open(args.infile, 'r') as f:
    for line in f:
#        print(line)
        line = line.rstrip()
        m_line = line[-22:]
        m_line = m_line[:-13]
        print(m_line)
        if not (m_line in bam_dict):
#            print(line)
            print(line, file=outfile)

