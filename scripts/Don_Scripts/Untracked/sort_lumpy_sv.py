#!/usr/bin/env python

from __future__ import print_function

import pytools
import argparse

parser = argparse.ArgumentParser(description='Sort Lumpy SV output')
parser.add_argument('infile', help='input lumpy output file to sort')
parser.add_argument('outfile', help='file to write sorted lumpy output to')

args = parser.parse_args()

outdata = []
chroms = { '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9':9, '10': 10, '11': 11, '12': 12,
           '13': 13, '14': 14, '15': 15, '16': 16, '17': 17, '18': 18, '19': 19, '20': 20, '21': 21, '22': 22,
           'X': 23, 'Y': 24, 'MT': 25 }

with open(args.infile) as f:
    for line in f:
#        print(str(outdata))
#        print(line)
        line = line.split()
        sv = pytools.parse_lumpy_line(line)
        if (not sv[13]['start_chrom'] in chroms) or (not sv[13]['end_chrom'] in chroms):
            continue
        for i in range(len(outdata)):
            if pytools.SV_before(sv, outdata[-i - 1]):
                continue
            else:
                if not i:
                    outdata.append(sv)
                    break
                else:
                    outdata.insert(-i, sv)
                    break
        else:
            outdata.insert(0, sv)

outfile = open(args.outfile, 'w')
for sv in outdata:
    evidence = ''
    if 'split' in sv[11]:
        if 'discordant' in sv[11]:
            evidence = '1,1;2,1'
        else:
            evidence = '2,1'
    else:
        evidence = '1,1'
    print('\t'.join(sv[:11] + ['IDS:' + evidence] + line[12:13] + ['MAX:' + sv[13]['start_chrom'] + ':' + str(sv[13]['start_pos']) + ';' + sv[13]['end_chrom'] + ':' + str(sv[13]['end_pos']) ] + [ '95:' + sv[14]['start_chrom'] + ':' + str(sv[14]['start_range'][0]) + '-' + str(sv[14]['start_range'][1]) + ';' + sv[14]['end_chrom'] + ':' + str(sv[14]['end_range'][0]) + '-' + str(sv[14]['end_range'][1]) ] ), file=outfile)
