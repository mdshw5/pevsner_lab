#!/usr/bin/env python

from __future__ import print_function

import argparse
import os
import subprocess, shlex

parser = argparse.ArgumentParser(description='Find valid control bam files and make a link to them in a second directory')
parser.add_argument('-inlist', help='A list of control files (one per line)')
parser.add_argument('-bam_directory', help='A directory to look for bam files in')
parser.add_argument('-link_directory', help='A directory to write soft links to')
args = parser.parse_args()

f = open('/dev/null', 'w')

with open(args.inlist) as f:
    for line in f:
        make_link = False
        print(line)
        #line format is .../SRR######.sra.ncbi_enc#
        bam = args.bam_directory + line[-23:-14] + '.bam'
        print(str(bam))
        #Check if file exists#
        if os.path.isfile(bam):
            print('Exists!')
            if os.stat(bam).st_size > 500000000:
                print('Is large!')
                command_line = "samtools idxstats " + bam
                shell_args = shlex.split(command_line)
                p = subprocess.Popen(shell_args, stdout=f, stderr=subprocess.PIPE)
                if not p.stderr.read():
                    print('EOF present!')
                    make_link = True
        #Make the link#
        if make_link:
            os.symlink(bam, args.link_directory + line[-23:-14] + '.bam')
