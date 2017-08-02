#!/usr/bin/env python3

import sys


case_file = open(sys.argv[1], 'r')

control_file = open(sys.argv[2], 'r')

#chroms = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'MT':25}

chroms = {'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24, 'chrM':25}

#sample_dict = {9:'01', 10:'02', 11:'03', 12:'04', 13:'05', 14:'01', 15:'02', 16:'03', 17:'04', 18:'05', 19:'06', 20:'06', 21:'07', 22:'07', 23:'08', 24:'08', 25:'09', 26:'09', 27:'10', 28:'10', 29:'11', 30:'11', 31:'12', 32:'12', 33:'13', 34:'13', 35:'14', 36:'14', 37:'15', 38:'15', 39:'16', 40:'16'}

options = []

if len(sys.argv) > 3:
    options = sys.argv[2:]

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
    case_line = case_file.readline()
    if case_line[0] == '#':
        print(case_line.rstrip())
        continue
    else:
        case_line = case_line.rstrip().split()
        break

while True:
    control_line = control_file.readline()
    if control_line[0] == '#':
        continue
    else:
        control_line = control_line.rstrip().split()
        break

while True:
    if not case_line[0] in chroms:
        case_line = case_file.readline().rstrip().split()
        if case_line == []:
            break
        continue

    elif not control_line[0] in chroms:
        control_line = control_file.readline().rstrip().split()
        if control_line == []:
            break
        continue

    elif (case_line[0] == control_line[0]) and (case_line[1] == control_line[1]):
        if control_line[6] == 'PASS':
            print('\t'.join(case_line))
        control_line = control_file.readline().rstrip().split()
        if control_line == []:
            break
        case_line = case_file.readline().rstrip().split()
        if case_line == []:
            break
        continue
    elif is_before(case_line, control_line):
        print('\t'.join(case_line))
        case_line = case_file.readline().rstrip().split()
        if case_line == []:
            break
        continue
    else:
        control_line = control_file.readline().rstrip().split()
        if control_line == []:
            break
        continue

