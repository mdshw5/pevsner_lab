#!/usr/bin/env python

import os
import subprocess
import shlex

input_files = os.listdir()

for input_file in input_files:
    base = input_file[:-4]
    outfile = open(base + '.sorted.gvf', 'w')
    command_line = 'srun /mnt/data/don/Scripts/gffsort.pl ' + input_file
    shell_args = shlex.split(command_line)
    p = subprocess.Popen(shell_args, stdout=subprocess.PIPE, universal_newlines=True)
    for line in p.stdout:
        line = line.rstrip()
        print(line, file=outfile)
    outfile.close()


for input_file in input_files:
    base = input_file[:-4]
    outfile = open(base + '.vat.gvf', 'w')
    command_line = 'srun -c 4 VAT --build hg19 -f /mnt/data/reference/VAAST/refGene_hg19_with_introns.header.sorted.gff3 -a /mnt/data/reference/hg19/ucsc.hg19.fasta ' + base + '.sorted.gvf'
    shell_args = shlex.split(command_line)
    p = subprocess.Popen(shell_args, stdout=subprocess.PIPE, universal_newlines=True)
    for line in p.stdout:
        line = line.rstrip()
        print(line, file=outfile)
    outfile.close()

command_line = 'srun VST -o "U(0..' + str(len(input_files) - 1) + ')" -b hg19 ' + ' '.join(input_files)
shell_args = shlex.split(command_line)
outfile = open('1KGenomes_control_exomes.cdr', 'w')
p = subprocess.Popen(shell_args, stdout=subprocess.PIPE, universal_newlines=True)
for line in p.stdout:
    line = line.rstrip()
    print(line, file=outfile)
