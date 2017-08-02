#!/usr/bin/env python3

import argparse
import os
import subprocess
import shlex

parser = argparse.ArgumentParser(description='Find which bams are incorrectly processed')
parser.add_argument('infiles', nargs='+', help='Input bam files to check if they are properly unpacked')
parser.add_argument('-outfile', help='A file to write incomplete bams to')
args = parser.parse_args()

outfile = open(args.outfile, 'w')
f = open('/dev/null', 'w')

for infile in args.infiles:
    replace = False
    if os.stat(infile).st_size < 500000000:
        replace = True
    if not replace:
        command_line = "samtools idxstats " + infile
        shell_args = shlex.split(command_line)
        p = subprocess.Popen(shell_args, stdout=f, stderr=subprocess.PIPE)
        for line in p.stderr:
            if len(line) > 10:
                replace = True
    if replace:
        print(str(infile), file=outfile)
