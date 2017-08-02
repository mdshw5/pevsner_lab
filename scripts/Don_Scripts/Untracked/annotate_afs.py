#!/usr/bin/env python

from __future__ import print_function

import sys

with open(sys.argv[1], 'r') as f:
    for line in f:
        line = line.rstrip()
        if line[0] == '#':
            print(line)
            continue
        difference = 0
        frequency = 0.5
        line = line.split()
        for genotype in line[9:]:
            genotype = genotype.split(':')
            if genotype[0] == '0/1':
                reads = genotype[1].split(',')
                if (int(reads[0]) + int(reads[1])) == 0:
                    continue
                alt_freq = float(reads[1]) / (float(reads[1]) + float(reads[0]))
                diff = abs(alt_freq - 0.5)
                if diff > difference:
                    difference = diff
                    frequency = alt_freq
        line[6] = str(frequency)
        print('\t'.join(line))
