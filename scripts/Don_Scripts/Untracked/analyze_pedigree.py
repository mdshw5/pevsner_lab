#!/usr/bin/env python3

import argparse
import shlex
import subprocess
import sys
import os.path
import time
import shutil
import os

# MUST BE RUN FROM $HOME #

'''
Required modules:
python3
python2.7
bwa
samtools
bedtools
samblaster
lumpy
VAAST
perl_dependencies
tabix_bgzip
'''

# The next step is to turn every task into a function to allow for abbreviated runs #

JDK = '/cm/shared/plab/apps/JDK/jdk1.8.0_11/bin/java'
GATK = '/cm/shared/plab/apps/GATK/3.3-0/GenomeAnalysisTK.jar'
REFERENCE = '/mnt/data/reference/hs37d5.fa'
EMIT_CONF = '10'
CALL_CONF = '20'
DBSNP = '/mnt/data/reference/GATK_b37/dbsnp_138.b37.vcf'
ANNOTATIONS = ['BaseQualityRankSumTest', 'RMSMappingQuality', 'TandemRepeatAnnotator', 'QualByDepth', 'MappingQualityRankSumTest', 'ReadPosRankSumTest', 'FisherStrand', 'Coverage', 'MappingQualityZero', 'SpanningDeletions', 'StrandOddsRatio']
FILTERS = '-BQRS 5 -FSBias 5 -mapq 10 -MQRS 5 -ReadPosRS 5 -QD 10'
SPECIAL = ['#', '$', '%', '*', '\\', '?', '[', ']', '\'', '"', ';', '&', '(', ')', '|', '^', '<', '>', '\n', '\t', ' ']
GFF = '/mnt/data/reference/VAAST/refGene_hg19_with_introns.header.sorted.gff3'
HG19 = '/mnt/data/reference/hg19/ucsc.hg19.fasta'
AF_VCFS = ['/mnt/data/reference/ExAC/ExAC.r0.2.sites.vep.vcf.gz', '/mnt/data/reference/GATK_b37/1000G_phase1.snps.high_confidence.b37.vcf.gz']
PRESENT_VCFS = ['/mnt/data/reference/GATK_b37/dbsnp_138.b37.vcf.gz']
BACKGROUND_CDR = '/mnt/data/reference/VAAST/1KG_refGene_Dec2011_CGDiv_NHLBI_NoCall.cdr'

parser = argparse.ArgumentParser(description='Takes a list of bam files as input, realigns them using bwa mem, and calls variants using the GATK.')
parser.add_argument('--cores', type=int, help='The number of cores to use for realignment.')
parser.add_argument('--jobs', type=int, help='The max number of jobs to submit simultaneously.')
parser.add_argument('--intermediate_directory', help='A directory to place intermediate files.')
parser.add_argument('--bam_directory', help='A directory to place final alignment files.')
parser.add_argument('--vcf_directory', help='A directory to place final variant calls.')
parser.add_argument('--results', help='A directory to place final resutls.')
parser.add_argument('--analysis_type', default='exome', help='Type of analysis, exome/genome [exome].')
parser.add_argument('fam_file', help='An input file in .fam format.')
parser.add_argument('inputs', nargs='+', help='Input files in .bam format.')
args = parser.parse_args()

if not os.path.isdir(args.intermediate_directory):
    sys.exit("ERROR: The intermedite_directory argument is mandatory.")
if not os.path.isdir(args.bam_directory):
    sys.exit("ERROR: The bam_directory argument is mandatory.")
if not os.path.isdir(args.vcf_directory):
    sys.exit("ERROR: The vcf_directory argument is mandatory.")
if not os.path.isdir(args.results):
    sys.exit("ERROR: The results argument is mandatory.")

if args.analysis_type != 'exome' and args.analysis_type != 'genome':
    sys.exit("ERROR: The '--analysis_type' is unrecognized.\nPossible options are genome or exome.")

if not args.intermediate_directory.endswith('/'):
    args.intermediate_directroy += '/'
if not args.bam_directory.endswith('/'):
    args.bam_directory += '/'
if not args.vcf_directory.endswith('/'):
    args.vcf_directory += '/'
if not args.results.endswith('/'):
    args.results += '/'

if not args.cores:
    args.cores = 1

if not args.jobs:
    args.jobs = 6

if not os.path.isfile(args.fam_file):
    sys.exit("ERROR: The supplied pedigree/fam file does not exist.")
for infile in args.inputs:
    if not os.path.isfile(infile):
        sys.exit("ERROR: The supplied file " + infile + " does not exist.")

#fix for multiple probands
probands = []
pedigrees = {}
with open(args.fam_file, 'r') as f:
    for line in f:
#        print(line.rstrip(), file=sys.stderr)
        if line[0] == '#':
            continue
        line = line.rstrip().split()
        if not line[0] in pedigrees:
            pedigrees[line[0]] = {}
            pedigrees[line[0]] = {'affected':[], 'unaffected':[], 'other':[]}
        pedigrees[line[0]][line[1]] = {}
        if line[2] != '0':
            pedigrees[line[0]][line[1]]['father'] = line[2]
        else:
            pedigrees[line[0]][line[1]]['father'] = None
        if line[3] != '0':
            pedigrees[line[0]][line[1]]['mother'] = line[3]
        else:
            pedigrees[line[0]][line[1]]['mother'] = None
        if line[4] == '1':
            pedigrees[line[0]][line[1]]['sex'] = 'male'
        elif line[4] == '2':
            pedigrees[line[0]][line[1]]['sex'] = 'female'
        else:
            pedigrees[line[0]][line[1]]['sex'] = None
        if line[5] == '1':
            pedigrees[line[0]][line[1]]['affected'] = True
            pedigrees[line[0]]['affected'].append(line[1])
        elif line[5] == '0':
            pedigrees[line[0]][line[1]]['affected'] = False
            pedigrees[line[0]]['unaffected'].append(line[1])
        else:
            pedigrees[line[0]][line[1]]['affected'] = None
            pedigrees[line[0]]['other'].append(line[1])
#        print(str(pedigrees) + '\n\n', file=sys.stderr)

print('INFO: Pedigrees: ' + str(pedigrees), file=sys.stderr)

for fam in pedigrees:
    probands.append(pedigrees[fam]['affected'][0])

samples = {}

infile_base = []
process_bin = []
command_list = []

for infile in args.inputs:
    if infile.endswith('.bam'):
        base = infile[infile.rfind('/') + 1 : -4]
    else:
        base = infile[infile.rfind('/') + 1:]
    infile_base.append(base)
#'''
    print('INFO: Processing ' + base + ' at: ' + time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
    command_line = 'srun -c ' + str(args.cores) + ' /mnt/data/don/Scripts/realign_bam.sh ' + infile + ' ' + base + ' ' + base + ' ' + args.intermediate_directory + ' ' + str(args.cores)
    command_list.append(command_line)
    shell_args = shlex.split(command_line)
    if len(process_bin) > int(args.jobs):
        #check if any processes have finished#
#        print('Checking process bin of len: ' + str(len(process_bin)), file=sys.stderr) 
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

hc_infiles = [''] * len(infile_base)

for i in range(len(infile_base)):
    hc_infiles[i] = infile_base[i] + '_sorted_mem.bam'
    shutil.move(args.intermediate_directory + hc_infiles[i], args.bam_directory)
    shutil.move(args.intermediate_directory + hc_infiles[i] + '.bai', args.bam_directory)
    hc_infiles[i] = args.bam_directory + hc_infiles[i]

#HaplotypeCaller

infiles = '-I ' + ' -I '.join(hc_infiles)
print('INFO: Running HaplotypeCaller on: ' + infiles + ' at: ' + time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
#memory is scaled with cores#
memory = int(int(args.cores) / 32.0 * 120.0)
memory = '-Xmx' + str(memory) + 'g'
output = args.intermediate_directory + probands[0]
annotations = '-A ' + ' -A '.join(ANNOTATIONS)
command_line = 'srun -c ' + str(args.cores) + ' ' + JDK + ' ' + memory + ' -jar ' + GATK + ' -T HaplotypeCaller -R ' + REFERENCE + ' -nct ' + str(args.cores) + ' --genotyping_mode DISCOVERY -stand_emit_conf ' + EMIT_CONF + ' -stand_call_conf ' + CALL_CONF + ' -D ' + DBSNP + ' -o ' + output + '_raw.vcf ' + annotations + ' ' + infiles
print('INFO: RUNNING: ' + command_line, file=sys.stderr)
shell_args = shlex.split(command_line)
p = subprocess.Popen(shell_args)
p.wait()

#'''

#Variant filtration#

print('INFO: Filtering variants at: ' + time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
output = args.intermediate_directory + probands[0]
depth = ''
if args.analysis_type == 'genome':
    depth = ' -depth ' + str(len(infile_base) * 7)
else:
    depth = ' -depth ' + str(len(infile_base) * 25)
filters = FILTERS + depth
command_line = 'srun -c 2 /mnt/data/don/Scripts/filter_vcf.py ' + filters + ' ' + output + '_raw.vcf ' + output + '_annotated.vcf'
print('INFO: RUNNING: ' + command_line, file=sys.stderr)
shell_args = shlex.split(command_line)
p = subprocess.Popen(shell_args)
p.wait()

command_line = "grep '^#\|PASS' " + output + '_annotated.vcf'
shell_args = shlex.split(command_line)
vcf_out = open(args.vcf_directory + probands[0] + '_pooled.vcf', 'w')
p = subprocess.Popen(shell_args, stdout=vcf_out)
p.wait()
vcf_out.close()

#bgzip/Tabix index#

print('INFO: Creating a bgzip/tabix indexed vcf at: ' + time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
zipped = open(args.vcf_directory + probands[0] + '_pooled.vcf.gz', 'w')
command_line = 'bgzip -c ' + args.vcf_directory + probands[0] + '_pooled.vcf'
shell_args = shlex.split(command_line)
print('INFO: RUNNING: ' + command_line, file=sys.stderr)
p = subprocess.Popen(shell_args, stdout=zipped)
p.wait()
zipped.close()

command_line = 'tabix -p vcf ' + args.vcf_directory + probands[0] + '_pooled.vcf.gz'
shell_args = shlex.split(command_line)
print('INFO: RUNNING: ' + command_line, file=sys.stderr)
p = subprocess.Popen(shell_args)
p.wait()

#VAAST#

print('INFO: Running vaast_converter at: ' + time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
command_line = 'srun -c 2 vaast_converter -n 0 --build hg19 --path ' + args.intermediate_directory + ' ' + args.vcf_directory + probands[0] + '_pooled.vcf'
print('INFO: RUNNING: ' + command_line, file=sys.stderr)
shell_args = shlex.split(command_line)
p = subprocess.Popen(shell_args)
p.wait()

#'''

gvf_bases = []
for proband in probands:
    if os.path.isfile(args.intermediate_directory + proband + '.gvf'):
        gvf_bases.append(args.intermediate_directory + proband + '.gvf')
    else:
        all_files = os.listdir(args.intermediate_directory)
        for possible in all_files:
            if possible.endswith('.gvf') and proband in possible:
                gvf_bases.append(possible)
                break
if len(gvf_bases) != len(probands):
    sys.exit('ERROR: Could not find gvf files for all probands')

process_bin = []
outfile_bin = []

for i, gvf_base in enumerate(gvf_bases):
    if '/' in gvf_base:
        gvf_bases[i] = gvf_base[gvf_base.rfind('/') + 1: ]
    gvf_bases[i] = gvf_bases[i][: -4]

#'''

for gvf_base in gvf_bases:
    sorted_gvf = open(args.intermediate_directory + gvf_base + '_sorted.gvf', 'w')
    print('INFO: Running gffsort at: ' + time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
    command_line = 'srun -c 2 /mnt/data/don/Scripts/gffsort.pl ' + args.intermediate_directory + gvf_base + '.gvf'
    print('INFO: RUNNING: ' + command_line, file=sys.stderr)
    shell_args = shlex.split(command_line)
    if len(process_bin) > int(args.jobs):
        #print('Checking process bin of len: ' + str(len(process_bin)), file=sys.stderr)
        while True:
            running = 0
            for process in process_bin:
                if process.poll() is None:
                    running += 1
            if running > int(args.jobs):
                time.sleep(300)
            else:
                break
    p = subprocess.Popen(shell_args, stdout=sorted_gvf)
    process_bin.append(p)
    outfile_bin.append(sorted_gvf)
for p in process_bin:
    p.wait()
for outfile in outfile_bin:
    outfile.close()

process_bin = []
outfile_bin = []

for gvf_base in gvf_bases:
    print('INFO: Running VAT at: ' + time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
    vat_file = open(args.intermediate_directory + gvf_base + '_vat.gvf', 'w')
    sex = ''
    for i in pedigrees:
        for j in pedigrees[i]['affected']:
            if j in gvf_base:
                sex = pedigrees[i][j]['sex']
    if not sex:
        sex = ''
    else:
        sex = ' --sex ' + sex
    command_line = 'srun -c 2 VAT --build hg19 -f ' + GFF + ' -a ' + HG19 + sex + ' ' + args.intermediate_directory + gvf_base + '_sorted.gvf'
    shell_args = shlex.split(command_line)
    if len(process_bin) > int(args.jobs):
        #print('Checking process bin of len: ' + str(len(process_bin)), file=sys.stderr)
        while True:
            running = 0
            for process in process_bin:
                if process.poll() is None:
                    running += 1
            if running > int(args.jobs):
                time.sleep(300)
            else:
                break
    print('INFO: RUNNING: ' + command_line, file=sys.stderr)
    p = subprocess.Popen(shell_args, stdout=vat_file)
    process_bin.append(p)
    outfile_bin.append(vat_file)
for p in process_bin:
    p.wait()
for outfile in outfile_bin:
    outfile.close()

#'''

process_bin = []
outfile_bin = []

for gvf_base in gvf_bases:
    print('INFO: Running VST at: ' + time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
    cdr_file = open(args.intermediate_directory + gvf_base + '.cdr', 'w')
    command_line = 'srun -c 2 VST -o "U(0..0)" -b hg19 ' + args.intermediate_directory + gvf_base + '_vat.gvf'
    shell_args = shlex.split(command_line)
    if len(process_bin) > int(args.jobs):
        #print('Checking process bin of len: ' + str(len(process_bin)), file=sys.stderr)
        while True:
            running = 0
            for process in process_bin:
                if process.poll() is None:
                    running += 1
            if running > int(args.jobs):
                time.sleep(300)
            else:
                break
    print('INFO: RUNNING: ' + command_line, file=sys.stderr)
    p = subprocess.Popen(shell_args, stdout=cdr_file)
    process_bin.append(p)
    outfile_bin.append(cdr_file)
for p in process_bin:
    p.wait()
for outfile in outfile_bin:
    outfile.close()

process_bin = []
outfile_bin = []

for gvf_base in gvf_bases:
    print('INFO: Running VAAST at: ' + time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
    command_line = 'srun -c ' + str(args.cores) + ' VAAST -gp 1e6 -iht n -m lrt -r 0.12 --use_aas_info y -e --indel --enable_splice_sites y -p ' + str(args.cores) + ' --no_max_allele_count -o ' + args.intermediate_directory + gvf_base + '.vaast ' + GFF + ' ' + BACKGROUND_CDR + ' ' + args.intermediate_directory + gvf_base + '.cdr'
    shell_args = shlex.split(command_line)
    if len(process_bin) > int(args.jobs):
        #print('Checking process bin of len: ' + str(len(process_bin)), file=sys.stderr)
        while True:
            running = 0
            for process in process_bin:
                if process.poll() is None:
                    running += 1
            if running > int(args.jobs):
                time.sleep(300)
            else:
                break
    print('INFO: RUNNING: ' + command_line, file=sys.stderr)
    p = subprocess.Popen(shell_args)
    process_bin.append(p)
for p in process_bin:
    p.wait()

process_bin = []
    
for gvf_base in gvf_bases:
    print('INFO: Annotating VAAST output at: ' + time.strftime("%Y-%m-%d %H:%M:%S"), file=sys.stderr)
    command_line = 'srun -c 1 /mnt/data/don/Scripts/Tracked/annotate_vaast.py ' + args.fam_file + ' ' + args.vcf_directory + probands[0] + '_pooled.vcf ' + args.intermediate_directory + gvf_base + '.vaast.simple ' + args.intermediate_directory + gvf_base + '_annotated.simple --af_vcf ' + ' '.join(AF_VCFS) + ' --presence ' + ' '.join(PRESENT_VCFS)
    shell_args = shlex.split(command_line)
    if len(process_bin) > int(args.jobs):
        while True:
            running = 0
            for process in process_bin:
                if process.poll() is None:
                    running += 1
            if running > int(args.jobs):
                time.sleep(300)
            else:
                break
    print('INFO: RUNNING: ' + command_line, file=sys.stderr)
    p = subprocess.Popen(shell_args)
    process_bin.append(p)
for p in process_bin:
    p.wait()


#process_bin = []
#out_array = []
 

# Running LUMPY is below
'''
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
'''
