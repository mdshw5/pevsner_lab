#!/usr/bin/env python3

import argparse
import subprocess
import shlex
import os
import sys
import math
import time

parser = argparse.ArgumentParser(description='Run the GATK HaplotypeCaller over multiple intervals and merge the output to a final file. Splits processing using the -L argument')
parser.add_argument("dispatch_cmd", help="The dispatch command to be used to run the job ['']")
parser.add_argument("cores", type=int, help="The number of cores to split the command across")
parser.add_argument("gatk", help="The gatk .jar file")
parser.add_argument("reference", help="The reference file to call variants with")
parser.add_argument("args", help="Arguments to the GATK")
parser.add_argument("output", help="The name of the output file to write to.")
parser.add_argument("--memory", help="Memory per job. [1g]")
args = parser.parse_args()

if not args.memory:
    args.memory = '1g'

if os.path.isfile(args.reference + '.fai'):
    fai = args.reference + '.fai'
else:
    sys.exit("No .fai file found for fasta file: " + args.reference + ".\nPlease create one using samtools faidx.")

total_size = 0
bp = []
chroms = []

with open(args.reference + '.fai', 'r') as f:
    for line in f:
        line = line.rstrip().split()
        total_size += int(line[1])
        chroms.append(line[0])
        bp.append(int(line[1]))

step_size = math.ceil(total_size / float(args.cores))
#print(str(step_size))

intervals = []
last_chrom = ''
last_pos = 1
passed = 0
cur_interval = []



for i in range(len(chroms)):
#    print("passed: " + str(passed))
#    print("pos: " + str(last_pos))
    if chroms[i] == 'hs37d5':
        continue
    while (last_pos + step_size - passed < bp[i]):
        cur_interval.append('-L ' + chroms[i] + ':' + str(last_pos) + '-' + str(last_pos + step_size - passed))
        intervals.append(' '.join(cur_interval))
        cur_interval = []
        last_pos += step_size - passed
        passed = 0
    else:
        cur_interval.append('-L ' + chroms[i] + ':' + str(last_pos) + '-' + str(bp[i]))
        passed += bp[i] - last_pos
        last_pos = 1

if cur_interval:
    intervals.append(' '.join(cur_interval))

'''
for interval in intervals:
    print(interval)
print(args.args)
'''

p = []
i = 0

if args.dispatch_cmd:
    args.dispatch_cmd += ' '

for interval in intervals:
    command = args.dispatch_cmd + 'java -Xmx' + args.memory + ' -Djava.io.tmpdir=/home/freedd/tmp  -jar ' + args.gatk + ' -T HaplotypeCaller -R ' + args.reference + ' -o tmp_' + str(i) + '.vcf ' + interval + ' ' + args.args
    print("Running: " + command)
    shell = shlex.split(command)
    p.append(subprocess.Popen(shell))
    time.sleep(20) # Can overload some job schedulers
    i += 1

variants = []

for i in range(len(intervals)):
    p[i].wait()
    variants.append('-V tmp_' + str(i) + '.vcf')

command = args.dispatch_cmd + 'java -Xmx' + args.memory + ' -cp ' + args.gatk + ' org.broadinstitute.gatk.tools.CatVariants -R ' + args.reference + ' ' + ' '.join(variants) + ' -out ' + args.output + ' -assumeSorted'
shell = shlex.split(command)
print("Running: " + command)
p = subprocess.Popen(shell)
p.wait()

for i in range(len(intervals)):
    os.remove('tmp_' + str(i) + '.vcf')
    os.remove('tmp_' + str(i) + '.vcf.idx')
