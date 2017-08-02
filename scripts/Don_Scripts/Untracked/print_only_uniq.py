#!/usr/bin/env python3

import sys

if len(sys.argv[1]) > 1:
    infile = open(sys.argv[1], 'r')
else:
    infile = sys.stdin

big_hash = {}

for line in infile:
    if line in big_hash:
        big_hash[line] = False
    else:
        big_hash[line] = True

for key, value in big_hash.items():
    if value:
        print(key)
