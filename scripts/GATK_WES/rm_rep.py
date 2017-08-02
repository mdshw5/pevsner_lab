#!/usr/bin/env python3

import argparse
import sys

def process_args():
    parser = argparse.ArgumentParser(description="Remove lines with a repeating first two columns")
    parser.add_argument("--infile", help="The input file")
    parser.add_argument("--outfile", help="The outut file")
    parser.add_argument("--skip_prefix", default='#', help="A line prefix to skip")
    return parser.parse_args()

def main(args):
    if not args:
        args = process_args()

    if args.infile:
        args.infile = open(args.infile, 'r')
    else:
        args.infile = sys.stdin

    if args.outfile:
        args.outfile = open(args.outfile, 'w')
    else:
        args.outfile = sys.stdout

    seen_before = set()

    for line in args.infile:
        line = line.rstrip()
        if line.startswith(args.skip_prefix):
            print(line, file=args.outfile)
            continue
        line = line.split()
        if line[0] + ':' + line[1] in seen_before:
            continue
        else:
            seen_before.add(line[0] + ':' + line[1])
            print('\t'.join(line), file=args.outfile)

if __name__ == "__main__":
    main(None)
