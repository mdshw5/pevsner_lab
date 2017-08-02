#!/usr/bin/env python3

import sys

for line in sys.stdin:
    line = line.rstrip().split()
    if len(line[3]) > 1:
        continue
    elif len(line[4]) > 1:
        continue
    else:
        print('\t'.join(line))
