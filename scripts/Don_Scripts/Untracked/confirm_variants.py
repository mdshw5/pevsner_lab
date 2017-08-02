#!/usr/bin/env python3

import argparse
import re

chroms = {
    '1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'MT':25, 'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24, 'chrM':25, 'GL000207.1':26, 'GL000226.1':27, 'GL000229.1':28, 'GL000231.1':29, 'GL000210.1':30, 'GL000239.1':31, 'GL000235.1':32, 'GL000201.1':33, 'GL000247.1':34, 'GL000245.1':35, 'GL000197.1':36, 'GL000203.1':37, 'GL000246.1':38, 'GL000249.1':39, 'GL000196.1':40, 'GL000248.1':41, 'GL000244.1':42, 'GL000238.1':43, 'GL000202.1':44, 'GL000234.1':45, 'GL000232.1':46, 'GL000206.1':47, 'GL000240.1':48, 'GL000236.1':49, 'GL000241.1':50, 'GL000243.1':51, 'GL000242.1':52, 'GL000230.1':53, 'GL000237.1':54, 'GL000233.1':55, 'GL000204.1':56, 'GL000198.1':57, 'GL000208.1':58, 'GL000191.1':59, 'GL000227.1':60, 'GL000228.1':61, 'GL000214.1':62, 'GL000221.1':63, 'GL000209.1':64, 'GL000218.1':65, 'GL000220.1':66, 'GL000213.1':67, 'GL000211.1':68, 'GL000199.1':69, 'GL000217.1':70, 'GL000216.1':71, 'GL000215.1':72, 'GL000205.1':73, 'GL000219.1':74, 'GL000224.1':75, 'GL000223.1':76, 'GL000195.1':77, 'GL000212.1':78, 'GL000222.1':79, 'GL000200.1':80, 'GL000193.1':81, 'GL000194.1':82, 'GL000225.1':83, 'GL000192.1':84
    }

parser = argparse.ArgumentParser(description='Using GIAB calls, evaluate a NA12878 variant callset')
parser.add_argument('GIAB_calls', help='A vcf file with variant calls from the GIAB consortium')
parser.add_argument('test_calls', help='Variant calls produced by the pipeline under evalutation')
parser.add_argument('bed', help='A bed file of regions to call variants over')
args = parser.parse_args()

giab_f = open(args.GIAB_calls, 'r')
test_f = open(args.test_calls, 'r')
bed = open(args.bed, 'r')

current_chrom = ''
current_start = ''
current_end = ''

# Initalize region #
line = bed.readline().rstrip().split()
current_chrom = line[0]
current_start = int(line[1])
current_end = int(line[2])

giab_line = ''
test_line = ''

# Initalize vcf lines #
while True:
    giab_line = giab_f.readline().rstrip()
    if giab_line[0] == '#':
        continue
    else:
        giab_line = giab_line.split()
        break

while True:
    test_line = test_f.readline().rstrip()
    if test_line[0] == '#':
        continue
    else:
        test_line = test_line.split()
        break

concordant = 0 # Called the same in GIAB and the pipeline
discordant = 0 # Called, but with a different zygosity ex. 0/1 to 1/1
not_called = 0 # Present in GIAB vcf, but not called in the test vcf
extra = 0      # Present in test vcf, but not called in the GIAB vcf

while True:
    if (chroms[current_chrom] < chroms[giab_line[0]] and chroms[current_chrom] < chroms[test_line[0]]): # Region is on an earlier chromosome
        line = bed.readline().rstrip().split()
        if not line:
            break
        current_chrom = line[0]
        current_start = int(line[1])
        current_end = int(line[2])
        continue
    if current_end < int(giab_line[1]) and current_end < int(test_line[1]): # End of region is at an earlier oposition
        line = bed.readline().rstrip().split()
        if not line:
            break
        current_chom = line[0]
        current_start = int(line[1])
        current_end = int(line[2])
        continue
    if chroms[test_line[0]] < chroms[current_chrom] or (chroms[test_line[0]] == chroms[current_chrom] and int(test_line[1]) < current_start): # Problematic with large indels
        test_line = test_f.readline().rstrip().split()
        if not test_line:
            break
        continue
    if chroms[giab_line[0]] < chroms[current_chrom] or (chroms[giab_line[0]] == chroms[current_chrom] and int(giab_line[1]) < current_start):# The giab vcf is behind. This may be problematic with large indels
        giab_line = giab_f.readline().rstrip().split()
        if not giab_line:
            break
        continue
    # Either one or both of the vcfs have a variant in the current region #
    if test_line[0] == giab_line[0] and test_line[1] == giab_line[1]: # Chromosome and position are the same
        giab_alleles = [giab_line[3]] + giab_line[4].split(',')
        test_alleles = [test_line[3]] + test_line[4].split(',')
        giab_geno = re.split('[/|]', giab_line[9].split(':')[0])
        test_geno = re.split('[/|]', test_line[9].split(':')[0])
        if (giab_alleles[int(giab_geno[0])] == test_alleles[int(test_geno[0])] and giab_alleles[int(giab_geno[1])] == test_alleles[int(test_geno[1])]) or \
                (giab_alleles[int(giab_geno[1])] == test_alleles[int(test_geno[0])] and giab_alleles[int(giab_geno[0])] == test_alleles[int(test_geno[1])]):
            print("Concordant: " + giab_line[0] + ':' + giab_line[1] + ' ' + giab_alleles[int(giab_geno[0])] + ' ' + giab_alleles[int(giab_geno[1])] + '\t' + test_line[0] + ':' + test_line[1] + ' ' + test_alleles[int(test_geno[0])] + ' ' + test_alleles[int(test_geno[1])])
            concordant += 1
            test_line = test_f.readline().rstrip().split()
            giab_line = giab_f.readline().rstrip().split()
            if not test_line or not giab_line:
                break
        elif (giab_alleles[int(giab_geno[0])] == test_alleles[int(test_geno[0])] or giab_alleles[int(giab_geno[0])] == test_alleles[int(test_geno[1])] or \
                giab_alleles[int(giab_geno[1])] == test_alleles[int(test_geno[0])] or giab_alleles[int(giab_geno[1])] == test_alleles[int(test_geno[1])]):
            print("Discordant: " + giab_line[0] + ':' + giab_line[1] + ' ' + giab_alleles[int(giab_geno[0])] + ' ' + giab_alleles[int(giab_geno[1])] + '\t' + test_line[0] + ':' + test_line[1] + ' ' + test_alleles[int(test_geno[0])] + ' ' + test_alleles[int(test_geno[1])])
            discordant += 1
            test_line = test_f.readline().rstrip().split()
            giab_line = giab_f.readline().rstrip().split()
            if not test_line or not giab_line:
                break
        else:
            print("Bad: " + giab_line[0] + ':' + giab_line[1] + ' ' + giab_alleles[int(giab_geno[0])] + ' ' + giab_alleles[int(giab_geno[1])] + '\t' + test_line[0] + ':' + test_line[1] + ' ' + test_alleles[int(test_geno[0])] + ' ' + test_alleles[int(test_geno[1])])
            not_called += 1
            extra += 1
            test_line = test_f.readline().rstrip().split()
            giab_line = giab_f.readline().rstrip().split()
            if not test_line or not giab_line:
                break
    else:
        if chroms[test_line[0]] < chroms[giab_line[0]] or (test_line[0] == giab_line[0] and int(test_line[1]) < int(giab_line[1])):
            test_alleles = [test_line[3]] + test_line[4].split(',')
            test_geno = re.split('[/|]', test_line[9].split(':')[0])
            extra += 1
            print("Extra: " + test_line[0] + ':' + test_line[1] + ' ' + test_alleles[int(test_geno[0])] + ' ' + test_alleles[int(test_geno[1])])
            test_line = test_f.readline().rstrip().split()
            if not test_line:
                break
        else:
            not_called += 1
            giab_alleles = [giab_line[3]] + giab_line[4].split(',')
            giab_geno = re.split('[/|]', giab_line[9].split(':')[0])
            print("NotCalled: " + giab_line[0] + ':' + giab_line[1] + ' ' + giab_alleles[int(giab_geno[0])] + ' ' + giab_alleles[int(giab_geno[1])])
            giab_line = giab_f.readline().rstrip().split()
            if not giab_line:
                break
print('\n\n')
print("Concordant: " + str(concordant))
print("Discordant: " + str(discordant))
print("Extra: " + str(extra))
print("Not called: " + str(not_called))
