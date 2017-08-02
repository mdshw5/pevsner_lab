#!/usr/bin/env python

from __future__ import print_function

import argparse

parser = argparse.ArgumentParser(description='Find sra runs with specific attributes. Print the files with those attributes to stdout.')
parser.add_argument('--attributes', help='A tab delimited sra attributes file')
parser.add_argument('--list', help='A list of SRR identifiers')
parser.add_argument('required', nargs='+', help='required attributes')
args = parser.parse_args()

big_dict = {}

with open(args.attributes, 'r') as f:
    for line in f:
        big_dict[line.rstrip().split()[8]] = line

with open(args.list, 'r') as f:
    for line in f:
        line = line.rstrip()
#        print(line)
        print_out = True
#        print(big_dict[line])
        for required in args.required:
#            print(attribute)
            if not required in big_dict[line]:
                print_out = False
#                print(big_dict[line])
        if print_out:
            print(line)
#            print(big_dict[line])
