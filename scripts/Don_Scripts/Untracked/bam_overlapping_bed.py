#!/usr/bin/env python3

import argparse
import subprocess
import shlex
import os.path

parser = argparse.ArgumentParser(description='From multiple bams, extract regions overlapping a bed file.')
parser.add_argument('bedfile', help='An input bed file of regions to extract')
parser.add_argument('bam_files', nargs='+', help='Bam files to extract reads from')
args = parser.parse_args()

with open(args.bedfile) as f:
    for line in f:
        line = line.rstrip().split()
        for bam in args.bam_files:
            outname = '/home/freedd/tmp/Patient_bams/' + os.path.basename(bam)[:-4] + '_' + line[0] + '_' + line[1] + '_' + line[2] + '.bam'
            command_line = 'samtools view -b -o ' + outname + ' ' + bam + ' ' + line[0] + ':' + line[1] + '-' + line[2]
            shell_args = shlex.split(command_line)
            p = subprocess.Popen(shell_args)
            p.wait()
