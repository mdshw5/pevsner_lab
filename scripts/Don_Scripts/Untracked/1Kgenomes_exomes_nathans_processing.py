#!/usr/bin/env python

from __future__ import print_function

import argparse
import subprocess
import shlex
import sys
import os.path
import re
import time
import os

parser = argparse.ArgumentParser(description='Process 1KGenomes Illumina Exome 90+bp paired-end reads through the GATK')
parser.add_argument('sequence_index', help='The sequence index from 1000 Genomes')
parser.add_argument('pedigree_file', help="The current complete pedigree information")
parser.add_argument('gvf_directory', help='A directory to store the GVF files')
args = parser.parse_args()

study_dict = {} #Individual ID + Study ID  #line[9] + line[3] #array:[0] = Individual ID [1] = array of read_ids
individual_dict = {} #Individual ID        #line[9]

ped_dict = {} #Individual ID   #line[9] #array: [0] = Array of parent ids; [1] = array of sib ids; [2] = Use (True/False)

pops = ['CEU', 'TSI', 'FIN', 'GBR', 'IBS']

with open(args.pedigree_file, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
        if not line[6] in pops:
            continue
        parents = []
        sibs = []
        use = True
        if line[2] != '0':
            parents.append(line[2])
        if line[3] != '0':
            parents.append(line[3])
        if line[8] != '0':
            if ',' in line[8]:
                sibs = line[8].translate(None, '"').split(',')
            else:
                sibs = [ line[8] ]
        ped_dict[line[1]] = [ parents, sibs, use ]

with open(args.sequence_index, 'r') as f:
    for line in f:
        line = line.rstrip().split('\t')
        if line[0] == 'FASTQ_FILE':
            continue
        if not line[10] in pops:
            continue
        if not line[12] == 'ILLUMINA':
            continue
        if line[18] == 'SINGLE':
            continue
        if line[20] == '1':
            continue
        if not line[25] == 'exome':
            continue
        if (not '_1' in line[0]) and (not '_2' in line[0]):
            continue
        read_length = int(line[24]) / int(line[23])
        if read_length < 90:
            continue
#        print('\t'.join(line))
        read_file = line[0].split('/')[-1]
        if not os.path.isfile('/mnt/data/1KGenomes_Exomes/' + read_file):
            print("WARNING: " + read_file + "does not exist. Retry Download", file=sys.stderr)
            continue
        print("INFO: Checking md5sum of /mnt/data/1KGenomes_Exomes/" + read_file, file=sys.stderr)
#        command_line = 'srun md5sum /mnt/data/1KGenomes_Exomes/' + read_file
#        shell_args = shlex.split(command_line)
#        p = subprocess.Popen(shell_args, stdout=subprocess.PIPE)
#        output = p.communicate()
#        md5sum = output[0].split()[0]
#        if md5sum != line[1]:
#            print("WARNING: " + read_file + " does not match its md5sum. Retry Download", file=sys.stderr)
#            continue
        if line[9] + line[3] in study_dict:
            study_dict[line[9] + line[3]][1].append(read_file)
        else:
            if line[9] in individual_dict:
                print("WARNING: " + line[9] + " is present in multiple studies", file=sys.stderr)
                individual_dict[line[9]].append(line[3])
                study_dict[line[9] + line[3]] = [ line[9], [ read_file ]]
            else:
                individual_dict[line[9]] = [ line[3] ]
                study_dict[line[9] + line[3]] = [ line[9], [ read_file ]]

#Make sure that the individual has no parents or sibs in the ped dict
for individual in study_dict:
    name = study_dict[individual][0]
    if not name in ped_dict:
        sys.exit("ERROR: The individual is not in the pedigree file: " + name)
    if not ped_dict[name][2]:
        continue
    if ped_dict[name][0]: #Parents
        for parent in ped_dict[name][0]:
            if parent in individual_dict:
                #Do not use if parental sample is available#
                print("INFO: " + name + " has at least one parental exome availble, ignoring file in processing", file=sys.stderr)
                ped_dict[name][2] = False
    if ped_dict[name][1]: #Sibs
        for sib in ped_dict[name][1]:
            #Do not use sibs. Mark them#
            if not sib in ped_dict:
                print("WARNING: " + sib + " is a sibling who is not found in the pedigree file", file=sys.stderr)
                continue
            ped_dict[sib][2] = False

'''
#sample_count = 0
process_bin = []
for individual in study_dict:
    name = study_dict[individual][0]
    if not ped_dict[name][2]:
        print("WARNING: Individual skipped due to pedigree information " + name, file=sys.stderr)
        continue
    read_count = 0
    for read in study_dict[individual][1]:
        if '_2' in read:
            continue
        pair = re.sub('_1', '_2', read)
        read_id = name + str(read_count)
        print("INFO: Running process_fqs.sh on " + read + " and " + pair + " with read group id " + read_id + " and sample name " + name, file=sys.stderr)
        command_line = 'srun -c 20 /mnt/data/don/Scripts/process_fqs.sh /mnt/data/1KGenomes_Exomes/' + read + ' /mnt/data/1KGenomes_Exomes/' + pair + ' ' + ' ' + read_id + ' ' + name
        shell_args = shlex.split(command_line)
        p = subprocess.Popen(shell_args)
        process_bin.append(p)
        time.sleep(10) #Don't overload the job handler#
        read_count += 1
for process in process_bin:
    process.wait()
'''

analysis_groups = [] #Groups may contain > 40 bam files, but only from 40 individuals
count = 0
group = []
for individual in study_dict:
    for bam_file in study_dict[individual][1]:
        if '_2' in bam_file:
            continue
        bam_file = '/mnt/data/1KGenomes_Exomes/' + re.sub('_1.filt.fastq.gz', '_sorted_mem.bam', bam_file)
        group.append(bam_file)
    count += 1
    if bool(count / 40):
        print("INFO: Analysis group contains: " + str(group), file=sys.stderr)
        analysis_groups.append(group)
        group = []
        count = 0
print("INFO: Analysis group contains: " + str(group), file=sys.stderr)
analysis_groups.append(group)

count = 0
process_bin = []
for group in analysis_groups:
    infiles = '-I ' + ' -I '.join(group)
    print("INFO: Running TargetCreator on: " + infiles, file=sys.stderr)
    command_line = 'srun -c 15 /cm/shared/plab/apps/JDK/jdk1.8.0_11/bin/java -Xmx70g -jar /cm/shared/plab/apps/GATK/3.2-2/GenomeAnalysisTK.jar -nt 15 -T RealignerTargetCreator -R /mnt/data/reference/hs37d5.fa -o Indelrealigner' + str(count) + '.intervals --known /mnt/data/reference/GATK_b37/Mills_and_1000G_gold_standard.indels.b37.vcf --known /mnt/data/reference/GATK_b37/1000G_phase1.indels.b37.vcf ' + infiles
    shell_args = shlex.split(command_line)
    p = subprocess.Popen(shell_args)
    process_bin.append(p)
    time.sleep(10) #Don't overload handler#
    count += 1
for process in process_bin:
    process.wait()

count = 0
process_bin = []
for group in analysis_groups:
    for bam_file in group:
        out_file = re.sub('_sorted_mem.bam', '_realigned.bam', bam_file)
        print("INFO: Running IndelRealigner on: " + bam_file, file=sys.stderr)
        command_line = 'srun -c 3 /cm/shared/plab/apps/JDK/jdk1.8.0_11/bin/java -Xmx5g -jar /cm/shared/plab/apps/GATK/3.2-2/GenomeAnalysisTK.jar -T IndelRealigner -R /mnt/data/reference/hs37d5.fa -targetIntervals Indelrealigner' + str(count) + '.intervals --knownAlleles /mnt/data/reference/GATK_b37/Mills_and_1000G_gold_standard.indels.b37.vcf --knownAlleles /mnt/data/reference/GATK_b37/1000G_phase1.indels.b37.vcf -I ' + bam_file + ' -o ' + out_file
        shell_args = shlex.split(command_line)
        p = subprocess.Popen(shell_args)
        process_bin.append(p)
        time.sleep(10) #Job handler#
    count += 1
for process in process_bin:
    process.wait()

count = 0
process_bin = []
for group in analysis_groups:
    infiles = '-I ' + ' -I '.join(group)
    infiles = re.sub('sorted_mem.bam', '_realigned.bam', infiles)
    print("INFO: Running HaplotypeCaller on: " + infiles, file=sys.stderr)
    command_line = 'srun -c 15 -R /mnt/data/reference/hs37d5.fa -nct 15 --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 20 -D /mnt/data/reference/GATK_b37/dbsnp_137.b37.vcf -o variants_out' + str(count) + '.vcf -A BaseQualityRankSumTest -A RMSMappingQuality -A TandemRepeatAnnotator -A QualByDepth -A MappingQualityRankSumTest -A ReadPosRankSumTest -A FisherStrand -A Coverage -A MappingQualityZero -A SpanningDeletions ' + infiles
    shell_args = shlex.split(command_line)
    p = subprocess.Popen(shell_args)
    process_bin.append(p)
    time.sleep(10)
    count += 1
for process in process_bin:
    process.wait()

count = 0
process_bin = []
for group in analysis_groups:
    print("INFO: Filtering Variants for variants_out" + str(count) + '.vcf', file=sys.stderr)
    command_line = 'srun -c 1 /mnt/data/don/Scripts/filter_vcf.py -BQRS 4 -depth 400 -FSBias 4 -mapq 10 -MQRS 4 -ReadPosRS 4 -QD 10 variants_out' + str(count) + '.vcf variants_annotated' + str(count) + '.vcf'
    shell_args = shlex.split(command_line)
    p = subprocess.Popen(shell_args)
    process_bin.append(p)
    time.sleep(10)
    count += 1
for process in process_bin:
    process.wait()

count = 0
process_bin = []
for group in analysis_groups:
    print("INFO: Creating filtered vcfs", file=sys.stderr)
    command_line = "srun -c 1 grep '^#\|PASS' variants_annotated" + str(count) + '.vcf'
    shell_args = shlex.split(command_line)
    p = subprocess.Popen(shell_args)
    process_bin.append(p)
    time.sleep(10)
    count += 1
for process in process_bin:
    process.wait()

count = 0
process_bin = []
for group in analysis_groups:
    print("INFO: Converting variants to GVF format", file=sys.stderr)
    command_line = 'srun vaast_converter -n 10 --build hg19 --path ' + str(args.gvf_directory) + ' variants_annotated' + str(count) + '.vcf'
    shell_args = shlex.split(command_line)
    p = subprocess.Popen(shell_args)
    time.sleep(10)
    count += 1
for process in process_bin:
    process.wait()

files = os.listdir(args.gvf_directory)
for infile in files:
    base = re.sub('.gvf', '', infile)
    print("INFO: Sorting GVF: " + infile, file=sys.stderr)
    outfile = open(base + '.sorted.gvf', 'w')
    command_line = 'srun /mnt/data/don/Scripts/gffsort.pl ' + infile
    shell_args = shlex.split(command_line)
    p = subprocess.Popen(shell_args, stdout=subprocess.PIPE)
    for line in p.stdout:
        line = line.rstrip()
        print(line, file=outfile)
    outfile.close()

for infile in files:
    base = re.sub('.gvf', '', infile)
    print("INFO: Annotating GVF: " + base, file=sys.stderr)
    outfile = open(base + '.vat.gvf', 'w')
    command_line = 'srun VAT --build hg19 -f /mnt/data/reference/VAAST/refGene_hg19_with_introns.header.sorted.gff3 -a /mnt/data/reference/hg19/ucsc.hg19.fasta ' + base + '.sorted.gvf'
    shell_args = shlex.split(command_line)
    p = subprocess.Popen(shell_args, stdout=subprocess.PIPE)
    for line in p.stdout:
        line = line.rstrip()
        print(line, file=outfile)
    outfile.close()

print("INFO: Creating CDR file", file=sys.stderr)
file_number = len(files) - 1
command_line = 'srun VST -o "U(0..' + str(file_number) + ')" -b hg19 ' + ' '.join(files)
shell_args = shlex.split(command_line)
outfile = open('1KGenomes_control_exomes.cdr', 'w')
p = subprocess.Popen(shell_args, stdout=subprocess.PIPE)
for line in p.stdout:
    line = line.rstrip()
    print(line, file=outfile)
