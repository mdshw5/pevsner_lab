#!/usr/bin/env python

from __future__ import print_function

import sys


vcf_file = open(sys.argv[1], 'r')
dbSnp_file = open(sys.argv[2], 'r')

chroms = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'MT':25}

sample_dict = {9:'01', 10:'02', 11:'03', 12:'04', 13:'05', 14:'01', 15:'02', 16:'03', 17:'04', 18:'05', 19:'06', 20:'06', 21:'07', 22:'07', 23:'08', 24:'08', 25:'09', 26:'09', 27:'10', 28:'10', 29:'11', 30:'11', 31:'12', 32:'12', 33:'13', 34:'13', 35:'14', 36:'14', 37:'15', 38:'15', 39:'16', 40:'16'}

options = []

if len(sys.argv) > 3:
    options = sys.argv[3:]

# possible options are:
# sample - only mutations in single samples
# snvs   - only snvs, exclude indels
# short  - only indels shorter than 4 bp
# single - only mutations in single individuals

def is_before(vcf1, vcf2):
#    chroms = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'MT':25}
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
#    print(vcf_line)
    if vcf_line[0] == '#':
        continue
    else:
        vcf_line = vcf_line.rstrip().split()
        break

while True:
    dbSnp_line = dbSnp_file.readline()
    if dbSnp_line[0] == '#':
        continue
    else:
        dbSnp_line = dbSnp_line.rstrip().split()
        break

total = 0
count = 0

histo = [0] * 201

while True:
#    print(str(vcf_line))
    debug = False
    if vcf_line[1] == '139688697':
#        print("Debug line")
        debug = True
    if vcf_line[6] != 'PASS': #Filter low qual variants
#        print("here1")
        vcf_line = vcf_file.readline().rstrip().split()
        if vcf_line == []:
            break
        continue
    elif not vcf_line[0] in chroms: #Only look at cononical chroms
#        print("here2")
        vcf_line = vcf_file.readline().rstrip().split()
        if vcf_line == []:
            break
        continue
    elif not dbSnp_line[0] in chroms: #Only cononical chroms
#        print("here3")
        dbSnp_line = dbSnp_file.readline().rstrip().split()
        if dbSnp_line == []:
            break
        continue
    elif 'snvs' in options:
#        print("Here4")
        if (len(vcf_line[3]) > 1 or len(vcf_line[4]) > 1):
#            print("Not an snv")
            vcf_line = vcf_file.readline().rstrip().split()
            if vcf_line == []:
                break
            continue
    elif 'short' in options:
#        print("here5")
        if(len(vcf_line[3]) > 3 or len(vcf_line[4]) > 3):
            vcf_line = vcf_file.readline().rstrip().split()
            if vcf_line == []:
                break
            continue
#    print("here6")
    if (vcf_line[0] == dbSnp_line[0] and vcf_line[1] == dbSnp_line[1]):
#        print("VCF line == dbsnp line")
        low_depth = True
        alt_allele = 1
        for tag in vcf_line[7].split(';'):
            if tag[0:2] == 'DP':
                if debug:
                    print("DP tag found")
                if int(tag[3:]) > (32 * 20):
                    low_depth = False
            if tag[0:2] == 'AC':
                if debug:
                    print("AC tag found")
                counts = tag[3:]
                counts = counts.split(',')
                if debug:
                    print("counts is: " + str(counts))
                    print("len counts is: " + str(len(counts)))
                if len(counts) > 1:
                    vcf_line = vcf_file.readline().rstrip().split()
                    if vcf_line == []:
                        break
                    dbSnp_line = dbSnp_file.readline().rstrip().split()
                    if dbSnp_line == []:
                        break
                    continue
                    
                    
                    '''
                        max_count = 0
                        max_pos = 0
                        for i in range(len(counts)):
                            if counts[i] > max_count:
                                max_count = counts[i]
                                max_pos = i
                        alt_allele = max_pos + 1
                        '''
                    
        if low_depth:
            vcf_line = vcf_file.readline().rstrip().split()
            if vcf_line == []:
                break
            dbSnp_line = dbSnp_file.readline().rstrip().split()
            if dbSnp_line == []:
                break
            continue
        else:
#            print("Depth looks good")
            if 'sample' in options:
                first_het = False
                second_het = False
                for genotype in vcf_line[9:]:
                    genotype = genotype.split(':')
                    if genotype[0] == '0/1':
                        if not first_het:
                            first_het = True
                        else:
                            second_het = True
                            break
                if first_het and second_het:
                    vcf_line = vcf_file.readline().rstrip().split()
                    if vcf_line == []:
                        break
                    dbSnp_line = dbSnp_file.readline().rstrip().split()
                    if dbSnp_line == []:
                        break
                    continue
            if 'single' in options:
                other = ''
                next = False
                for i in range(len(vcf_line[9:])):
                    i += 9
                    genotype = vcf_line[i].split(':')
                    if genotype[0] == '0/1':
                        if other:
                            if not other == sample_dict[i]:
                                next = True
                        else:
                            other = sample_dict[i]
                if next:
                    vcf_line = vcf_file.readline().rstrip().split()
                    if vcf_line == []:
                        break
                    dbSnp_line = dbSnp_file.readline().rstrip().split()
                    if dbSnp_line == []:
                        break
                    continue
#            print("Examining genotypes")
            for genotype in vcf_line[9:]:
                genotype = genotype.split(':')
                if genotype[0] == '0/1':
                    allele_counts = genotype[1].split(',')
                    if len(allele_counts) == 1:
                        continue
                    ref_reads = allele_counts[0]
                    alt_reads = allele_counts[alt_allele]
                    if (int(ref_reads) + int(alt_reads)) > 20:
                        count += 1
                        result = int(alt_reads) / (float(ref_reads) + int(alt_reads)) * 100
                        if result < 2 or result > 98:
                            print(vcf_line)
                            print(genotype)
                        total += result
                        bin = 0
                        if ( (result % 1) > 0.5 ):
                            bin = -1
                        histo[int(round(result) * 2) + bin] += 1
#            print("Done reading genotypes")
            vcf_line = vcf_file.readline().rstrip().split()
            if vcf_line == []:
                break
            dbSnp_line = dbSnp_file.readline().rstrip().split()
            if dbSnp_line == []:
                break
            continue
            
    elif is_before(vcf_line, dbSnp_line):
#        print("Vcf line is before")
        vcf_line = vcf_file.readline().rstrip().split()
        if vcf_line == []:
            break
        continue
    else:
#        print("dbSnp line is before")
        dbSnp_line = dbSnp_file.readline().rstrip().split()
        if dbSnp_line == []:
            break
        continue
        
average = total / count
print("Count is: " + str(count))
print("Average is: " + str(average))

print(str(histo))
