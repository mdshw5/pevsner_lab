#!/usr/bin/env python3

import argparse
import sys
from pyfaidx import Fasta

'''
TO DO:
-Handle indels
-Handle splic sites
-Improve algorithm when switching chromosomes
'''

parser = argparse.ArgumentParser(description='Annotate variants that are in exons (including insertions/deletions)')
parser.add_argument('ref_genes', help='A file containing gene annotations')
parser.add_argument('ref_seq', help='A reference genome sequence')
parser.add_argument('--infile', help='An input file (default is stdin)')
args = parser.parse_args()

infile = None
if args.infile:
    infile = open(args.infile, 'r')
else:
    infile = sys.stdin

reference = Fasta(args.ref_seq)

chroms = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24, 'MT':25, 'chr1':1, 'chr2':2, 'chr3':3, 'chr4':4, 'chr5':5, 'chr6':6, 'chr7':7, 'chr8':8, 'chr9':9, 'chr10':10, 'chr11':11, 'chr12':12, 'chr13':13, 'chr14':14, 'chr15':15, 'chr16':16, 'chr17':17, 'chr18':18, 'chr19':19, 'chr20':20, 'chr21':21, 'chr22':22, 'chrX':23, 'chrY':24, 'chrM':25}

ref_file = open(args.ref_genes, 'r')
ref_line = ref_file.readline().rstrip().split()

codon = {     'ttt': 'F', 'tct': 'S', 'tat': 'Y', 'tgt': 'C',
              'ttc': 'F', 'tcc': 'S', 'tac': 'Y', 'tgc': 'C',
              'tta': 'L', 'tca': 'S', 'taa': '*', 'tga': '*',
              'ttg': 'L', 'tcg': 'S', 'tag': '*', 'tgg': 'W',
              'ctt': 'L', 'cct': 'P', 'cat': 'H', 'cgt': 'R',
              'ctc': 'L', 'ccc': 'P', 'cac': 'H', 'cgc': 'R',
              'cta': 'L', 'cca': 'P', 'caa': 'Q', 'cga': 'R',
              'ctg': 'L', 'ccg': 'P', 'cag': 'Q', 'cgg': 'R',
              'att': 'I', 'act': 'T', 'aat': 'N', 'agt': 'S',
              'atc': 'I', 'acc': 'T', 'aac': 'N', 'agc': 'S',
              'ata': 'I', 'aca': 'T', 'aaa': 'K', 'aga': 'R',
              'atg': 'M', 'acg': 'T', 'aag': 'K', 'agg': 'R',
              'gtt': 'V', 'gct': 'A', 'gat': 'D', 'ggt': 'G',
              'gtc': 'V', 'gcc': 'A', 'gac': 'D', 'ggc': 'G',
              'gta': 'V', 'gca': 'A', 'gaa': 'E', 'gga': 'G',
              'gtg': 'V', 'gcg': 'A', 'gag': 'E', 'ggg': 'G'
        }
         

def reverse_compliment(seq):
    out_seq = ''
    for i in seq:
        if i == 'A':
            out_seq = 'T' + out_seq
        elif i == 'C':
            out_seq = 'G' + out_seq
        elif i == 'G':
            out_seq = 'C' + out_seq
        else:
            out_seq = 'A' + out_seq
    return out_seq

for line in infile:
    line = line.rstrip().split()
    before = False #Guarenteed to be before the variant#
    if line[0] == 'MT' or line[0] == 'chrM':
        #Skip mitochondrial variants#
        continue
    while True:
        if ref_line == []:
            ref_file.seek(0)
            ref_line = ref_file.readline().rstrip().split()
            before = True
            continue
        try:
            chroms[line[0]] == chroms[ref_line[2]]
        except IndexError:
            print(line, file=sys.stderr)
            print(ref_line, file=sys.stderr)
            raise
        if chroms[line[0]] == chroms[ref_line[2]]:
            if int(line[1]) < int(ref_line[4]):
                if before:
                    print('\t'.join(line))
                    break
                else:
                    if ref_file.tell() < 2001:
                        ref_file.seek(0)
                        ref_line = ref_file.readline().rstrip().split()
                    else:
                        ref_file.seek(ref_file.tell() - 2000, 0)
                        ref_line = ref_file.readline().rstrip().split()
                        ref_line = ref_file.readline().rstrip().split()
                    before = True
            else:
                if int(line[1]) <= int(ref_line[5]):
                    #Not intergenic; intronic or exonic#
                    if int(line[1]) <= int(ref_line[6]):
                        #5'UTR
                        if '5UTR=' + ref_line[12] in line[7]:
                            pass
                        else:
                            line[7] = '5UTR=' + ref_line[12] + ';' + line[7]
                        ref_line = ref_file.readline().rstrip().split()
                        continue
                    else:
                        if int(line[1]) >= int(ref_line[7]):
                            #3'UTR
                            if '3UTR=' + ref_line[12] in line[7]:
                                pass
                            else:
                                line[7] = '3UTR=' + ref_line[12] + ';' + line[7]
                            ref_line = ref_file.readline().rstrip().split()
                            continue
                        else:
                            #Intronic or coding#
                            coding = False
                            exon_start = ref_line[9].split(',')
                            exon_end = ref_line[10].split(',')
                            offset = ref_line[15].split(',')
                            for i in range(len(exon_start) - 1):
                                if int(line[1]) >= int(exon_start[i]) and int(line[1]) <= int(exon_end[i]):
                                    coding = True
                                    if offset[i] == '-1':
                                        if 'NONCODING=' + ref_line[12] in line[7]:
                                            pass
                                        else:
                                            line[7] = 'NONCODING=' + ref_line[12] + ';' + line[7]
                                        break
                                    else:
                                        if len(line[3]) > 1 or len(line[4]) > 1:
                                            #insertion or deletion#
                                            if len(line[3]) > 1:
                                                #deletion#
                                                if (len(line[3]) - 1) % 3 == 0:
                                                    #Non frameshift#
                                                    if 'NON_FRAME_DEL=' + ref_line[12] in line[7]:
                                                        pass
                                                    else:
                                                        line[7] = 'NON_FRAME_DEL=' + ref_line[12] + ';' + line[7]
                                                else:
                                                    #Frameshift#
                                                    if 'FRAMESHIFT_DEL=' + ref_line[12] in line[7]:
                                                        pass
                                                    else:
                                                        line[7] = 'FRAMESHIFT_DEL=' + ref_line[12] + ';' + line[7]
                                            else:
                                                #Insertion#
                                                if (len(line[4]) - 1) % 3 == 0:
                                                    #Non frameshift#
                                                    if 'NON_FRAME_INS=' + ref_line[12] in line[7]:
                                                        pass
                                                    else:
                                                        line[7] = 'NON_FRAME_INS=' + ref_line[12] + ';' + line[7]
                                                else:
                                                    #Frameshift#
                                                    if 'FRAMESHIFT_INS=' + ref_line[12] in line[7]:
                                                        pass
                                                    else:
                                                        line[7] = 'FRAMESHIFT_INS=' + ref_line[12] + ';' + line[7]
                                                    
                                        else:
                                            if ref_line[3] == '+':
                                                adjust = ( (int(line[1]) - int(exon_start[i]) - 1) - (3 - int(offset[i])) ) % 3
                                                start = int(line[1]) - adjust - 1
                                                end = start + 3
                                                ref_seq = reference[line[0]][start:end].seq
                                                if adjust == 0:
                                                    var_seq = line[4] + ref_seq[1:]
                                                elif adjust == 1:
                                                    var_seq = ref_seq[0] + line[4] + ref_seq[2]
                                                elif adjust == 2:
                                                    var_seq = ref_seq[:2] + line[4]
                                                else:
                                                    sys.exit("Adjust not between 0 and 2.")
                                            else:
                                                adjust = (int(exon_end[i]) + int(offset[i]) - int(line[1])) % 3
                                                start = int(line[1]) - (2 - adjust) - 1
                                                end = start + 3
                                                ref_seq = reference[line[0]][start:end].seq
                                                if adjust == 0:
                                                    var_seq = ref_seq[:2] + line[4] 
                                                elif adjust == 1:
                                                    var_seq = ref_seq[0] + line[4] + ref_seq[2]
                                                elif adjust == 2:
                                                    var_seq = line[4] + ref_seq[1:]
                                                else:
                                                    sys.exit("Adjust not between 0 and 2.")
                                                ref_seq = reverse_compliment(ref_seq)
                                                var_seq = reverse_compliment(var_seq)
                                            '''
                                            adjust = abs( (int(exon_end[i]) - int(line[1]) - 2) - (3 - int(offset[i])) + 1 ) % 3
                                            start = int(line[1]) - adjust - 1
                                            end = start + 3
                                            ref_seq = reference[line[0]][start:end].seq
                                            if adjust == 0:
                                                var_seq = line[4] + ref_seq[1:]
                                            elif adjust == 1:
                                                var_seq = ref_seq[0] + line[4] + ref_seq[2]
                                            elif adjust == 2:
                                                var_seq = ref_seq[:2] + line[4]
                                            else:
                                                sys.exit("Adjust not between 0 and 2.")
                                            ref_seq = reverse_compliment(ref_seq)
                                            var_seq = reverse_compliment(var_seq)
                                            '''
                                            
                                            if codon[ref_seq.lower()] == codon[var_seq.lower()]:
                                                if 'SYN=' + ref_line[12] in line[7]:
                                                    pass
                                                else:
                                                    line[7] = 'SYN=' + ref_line[12] + ';' + line[7]
                                            elif codon[var_seq.lower()] == '*':
                                                if 'NONSENSE=' + ref_line[12] in line[7]:
                                                    pass
                                                else:
                                                    line[7] = 'NONSENSE=' + ref_line[12] + ';' + line[7]
                                            else:
                                                if 'MISSENSE=' + codon[ref_seq.lower()] + '->' + codon[var_seq.lower()] + ',' + ref_line[12] in line[7]:
                                                    pass
                                                else:
                                                    line[7] = 'MISSENSE=' + codon[ref_seq.lower()] + '->' + codon[var_seq.lower()] + ',' + ref_line[12] + ';' + line[7]
                                            break
                            if not coding:
                                if 'INTRONIC=' + ref_line[12] in line[7]:
                                    pass
                                else:
                                    line[7] = 'INTRONIC=' + ref_line[12] + ';' + line[7]
                            ref_line = ref_file.readline().rstrip().split()
                            continue
                else:
                    before = True
                    ref_line = ref_file.readline().rstrip().split()
                    continue
        else:
            if before:
                print('\t'.join(line))
                break
            ref_line = ref_file.readline().rstrip().split()
        
    
