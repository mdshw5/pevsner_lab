#!/usr/bin/env python

from __future__ import print_function

import sys


vcf_file = open(sys.argv[1], 'r')

chroms = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'MT':25}

sample_dict = {9:'01', 10:'02', 11:'03', 12:'04', 13:'05', 14:'01', 15:'02', 16:'03', 17:'04', 18:'05', 19:'06', 20:'06', 21:'07', 22:'07', 23:'08', 24:'08', 25:'09', 26:'09', 27:'10', 28:'10', 29:'11', 30:'11', 31:'12', 32:'12', 33:'13', 34:'13', 35:'14', 36:'14', 37:'15', 38:'15', 39:'16', 40:'16'}

options = []

if len(sys.argv) > 2:
    options = sys.argv[2:]

# snvs   - only mutations which are SNVs
# short  - only SNVs and indels 3bp or less
# sample - only mutations in a single sample
# single - only mutations in single individuals
# AFs    - only mutations with allele frequencies between 0 and 0.3 and 0.7 and 1

def is_before(vcf1, vcf2):
    if chroms[vcf1[0]] < chroms[vcf2[0]]:
        return True
    elif chroms[vcf1[0]] > chroms[vcf2[0]]:
        return False
    else:
        if int(vcf1[1]) <= int(vcf2[1]):
            return True
        else:
            return False

#Remove headers
while True:
    vcf_line = vcf_file.readline()
    if vcf_line[0] == '#':
        print(vcf_line.rstrip())
        continue
    else:
        vcf_line = vcf_line.rstrip().split()
        break


while True:
    if vcf_line[6] != 'PASS': #Filter low qual variants
        vcf_line = vcf_file.readline().rstrip().split()
        if vcf_line == []:
            break
        continue
    elif not vcf_line[0] in chroms: #Only look at cononical chroms
        vcf_line = vcf_file.readline().rstrip().split()
        if vcf_line == []:
            break
        continue
    elif 'snvs' in options:
        if (len(vcf_line[3]) > 1 or len(vcf_line[4]) > 1):
            vcf_line = vcf_file.readline().rstrip().split()
            if vcf_line == []:
                break
            continue
    elif 'short' in options:
        if(len(vcf_line[3]) > 3 or len(vcf_line[4]) > 3):
            vcf_line = vcf_file.readline().rstrip().split()
            if vcf_line == []:
                break
            continue
    next_line = False
    alt_allele = 1
    for tag in vcf_line[7].split(';'):
        if tag[0:2] == 'DP':
            if int(tag[3:]) < (32 * 20):
                next_line = True
        if tag[0:2] == 'AC':
            counts = tag[3:]
            counts = counts.split(',')
            if len(counts) > 1:
                next_line = True
                    
    if next_line:
        vcf_line = vcf_file.readline().rstrip().split()
        if vcf_line == []:
            break
        continue
    if 'sample' in options:
        first_het = False
        second_het = False
        for genotype in vcf_line[9:]:
            genotype = genotype.split(':')
            if genotype[0] == '0/1' or genotype[0] == '1/1':
                if not first_het:
                    first_het = True
                else:
                    second_het = True
                    break
        if first_het and second_het:
            vcf_line = vcf_file.readline().rstrip().split()
            if vcf_line == []:
                break
            continue
    if 'single' in options:
        other = ''
        next = False
        for i in range(len(vcf_line[9:])):
            i += 9
            genotype = vcf_line[i].split(':')
            if genotype[0] == '0/1' or genotype[0] == '1/1':
                if other:
                    if not other == sample_dict[i]:
                        next = True
                else:
                    other = sample_dict[i]
        if next:
            vcf_line = vcf_file.readline().rstrip().split()
            if vcf_line == []:
                break
            continue
    if 'sample' in options:
        print('\t'.join(vcf_line))
    if 'AFs' in options:
        print_line = False
        for i in range(len(vcf_line[9:])):
            i += 9
            genotype = vcf_line[i].split(':')
            if genotype[0] == '0/1':
                allele_counts = genotype[1].split(',')
                if len(allele_counts) == 1:
                    continue
                ref_reads = allele_counts[0]
                alt_reads = allele_counts[alt_allele]
                if (int(ref_reads) + int(alt_reads)) > 20:
                    result = int(alt_reads) / (float(ref_reads) + int(alt_reads)) * 100
                    if result <= 30 or result >= 70:
                        print_line = True
        if print_line:
            print('\t'.join(vcf_line))
    vcf_line = vcf_file.readline().rstrip().split()
    if vcf_line == []:
        break
    continue
