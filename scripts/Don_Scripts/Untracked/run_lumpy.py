#!/usr/bin/env python3

import argparse
import shlex
import subprocess
import sys
import os.path
import time

parser = argparse.ArgumentParser(description='Takes a list of bam files as input, realigns them using bwa mem, and calls variants using LUMPY.')
parser.add_argument('--filter_file', help='A file used to define which inputs to process')
parser.add_argument('--filter_sample_column', help='In the filter file, a column which as the sample id. Sample id must = input file base name.')
parser.add_argument('--filter_text', nargs='+', help='Text which must be found in the sample line of the filter file for that sample to be output.')
parser.add_argument('--cores', type=int, help='The number of cores to use for realignment.')
parser.add_argument('--jobs', type=int, help='The max number of jobs to submit simultaneously.')
parser.add_argument('intermediate_directory', help='A directory to place intermediate files.')
parser.add_argument('output_directory', help='A directory to output .lumpy_out.txt files to.')
parser.add_argument('inputs', nargs='+', help='Input files in .bam format.')
args = parser.parse_args()

if not os.path.isdir(args.intermediate_directory):
    sys.exit("The first positional argument must be a directory")
if not os.path.isdir(args.output_directory):
    sys.exit("The second positional argument must be a directory")

if not args.intermediate_directory.endswith('/'):
    args.intermediate_directroy += '/'
if not args.output_directory.endswith('/'):
    args.output_directory += '/'

if not args.cores:
    sys.exit('Cores is a required argument')

if not args.jobs:
    sys.exit('Jobs is a required argument')

samples = {}

if args.filter_file:
    if not args.filter_sample_column or not args.filter_text:
        sys.exit("If a filtering file is supplied --filter_sample_column and --filter_text are required arguments")
    with open(args.filter_file, 'r') as f:
        for line in f:
            print_line = [False] * len(args.filter_text)
            line = line.rstrip().split('\t')
            for column in line:
                for i in range(len(args.filter_text)):
                    if column == args.filter_text[i]:
                        print_line[i] = True
            samples[line[int(args.filter_sample_column)]] = all(print_line)

infile_base = []
process_bin = []
command_list = []

for infile in args.inputs:
    if infile.endswith('.bam'):
        base = infile[infile.rfind('/') + 1:-4]
    else:
        print('Input file does not end with ".bam". Is the input file in bam format? Skipping: ' + infile, file=sys.stderr)
        continue
    if args.filter_file:
        if base in samples:
            if samples[base]:
                infile_base.append(base)
            else:
                print(base + ' does not meet filter. Skipping.', file=sys.stderr)
                continue
        else:
            print(base[-1] + ' is not found in filter file. Skipping.', file=sys.stderr)
            continue
    else:
        infile_base.append(base)
    if os.path.isfile(args.intermediate_directory + base + '_sr_sorted.bam'):
        print('Bam is already realigned. Skipping: ' + infile_base[-1], file=sys.stderr)
        continue
    else:
        print('Processing ' + infile_base[-1], file=sys.stderr)
    command_line = 'srun -c ' + str(args.cores) + ' /mnt/data/don/Scripts/realign_bam.sh ' + infile + ' ' + infile_base[-1] + ' ' + infile_base[-1] + ' ' + args.intermediate_directory + ' ' + str(args.cores)
    command_list.append(command_line)
    shell_args = shlex.split(command_line)
    if len(process_bin) > int(args.jobs):
        #check if any processes have finished#
        print('Checking process bin of len: ' + str(len(process_bin)), file=sys.stderr) 
        while True:
            running = 0
            for process in process_bin:
                if process.poll() is None:
                    running += 1
            if running > int(args.jobs):
                time.sleep(300)
            else:
                break
    p = subprocess.Popen(shell_args)
    process_bin.append(p)
    time.sleep(60)
for process in process_bin:
    process.wait()


process_bin = []
out_array = []

for base in infile_base:
    #Get read length#
    read_length = ''
    command_line = 'samtools view ' + args.intermediate_directory + base + '_sorted_mem.bam'
    shell_args = shlex.split(command_line)
    p = subprocess.Popen(shell_args, stdout=subprocess.PIPE, universal_newlines=True)
    for line in p.stdout:
        line = line.rstrip().split()
        if 'H' in line[5]:
            continue
        else:
            read_length = str(len(line[9]))
            p.terminate()
            break
    #Get the paired end distribution info#
    distro_file = open(args.intermediate_directory + base + '.distro', 'r')
    distro = distro_file.readline().rstrip().split()
    mean = distro[0][5:]
    stdev = distro[1][6:]
    distro_file.close()
    command_line = 'srun -c 4 lumpy -mw 3 -tt 0.0 -sr bam_file:' + args.intermediate_directory + base + '_sr_sorted.bam,back_distance:40,min_mapping_threshold:12,weight:1,id:1 -pe bam_file:' + args.intermediate_directory + base + '_disc_sorted.bam,histo_file:' + args.intermediate_directory + base + '.pe.histo,mean:' + mean + ',stdev:' + stdev + ',read_length:' + read_length + ',min_non_overlap:' + str(int(read_length) - 10) + ',discordant_z:2,back_distance:20,min_mapping_threshold:12,weight:1,id:2'
    shell_args = shlex.split(command_line)
    out_array.append(open(args.output_directory + base + '.lumpy_out.txt', 'w'))
    if len(process_bin) > int(args.jobs):
        #check if any processes have finished#
        print('Checking process bin of len: ' + str(len(process_bin)), file=sys.stderr)
        while True:
            running = 0
            for i in range(len(process_bin)):
                if process_bin[i].poll() is None:
                    running += 1
            if running > int(args.jobs):
                time.sleep(300)
            else:
                break
    p = subprocess.Popen(shell_args, stdout=out_array[-1])
    process_bin.append(p)
    time.sleep(30)

for process in process_bin:
    process.wait()
