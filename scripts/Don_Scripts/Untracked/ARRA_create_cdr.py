#!/usr/bin/env python

from __future__ import print_function

import argparse
import subprocess
import shlex
import os
import sys

parser = argparse.ArgumentParser(description='In a directory of GVF files, process only the files which have autism')
parser.add_argument('gvf_directory', help='A directory containing *.vat.gvf files')
parser.add_argument('status_file', help='A file with subject id in the second column, consent status in the third and affected status in the fourth.')
parser.add_argument('translate', help='A file to translate sampleID (column 3) to subjectID (column 1)')
parser.add_argument('outfile', help='A file to write output to')
args = parser.parse_args()

translate_dict = {}

ped_dict = {} #Key = subject ID. Value is a (consent, affected_status)

with open(args.translate, 'r') as f:
    for line in f:
        if line[0] == '#':
            continue
        line = line.rstrip().split('\t')
        translate_dict[line[2]] = line[0]

with open(args.status_file, 'r') as f:
    for line in f:
        if line == [''] or line == '\t' or line == '\n':
            continue
        if line[0] == '#':
            continue
        line = line.rstrip().split('\t')
        if line[0] == 'dbGaP SubjID':
            continue
        try:
            ped_dict[line[1]] = ( line[2], line[3] )
        except IndexError:
            print(str(line), file=sys.stderr)
            sys.stderr.flush()
            sys.exit("IndexError")

aut_files = os.listdir(args.gvf_directory)

run_files = []

for aut_file in aut_files:
    if aut_file[-8:] != '.vat.gvf':
        continue
    sample = aut_file[:-8]
    if sample in translate_dict:
        sample = translate_dict[sample]
    else:
        print("WARNING: " + sample + " is not in the translate dictionary", file=sys.stderr)
        continue
    if sample in ped_dict:
        if ped_dict[sample][0] == '2':
            print("WARNING: " + sample + " has consent for GRU", file=sys.stderr)
        elif ped_dict[sample][1] == '1':
            print("WARNING: " + sample + " is unaffected and is not related to an ASD individual", file=sys.stderr)
        elif ped_dict[sample][1] == '3':
            print("INFO: " + sample + " is an unaffected family member", file=sys.stderr)
        elif ped_dict[sample][1] == '2':
            run_files.append(args.gvf_directory + '/' + aut_file)
    else:
        print("WARNING: " + sample + " has no representation in the pedigree dictionary", file=sys.stderr)
        continue

file_number = len(run_files) - 1
command_line = 'srun VST -o "U(0..' + str(file_number) + ')" -b hg19 ' + ' '.join(run_files)
shell_args = shlex.split(command_line)
outfile = open(args.outfile, 'w')
p = subprocess.Popen(shell_args, stdout=subprocess.PIPE)
for line in p.stdout:
    line = line.rstrip()
    print(line, file=outfile)
