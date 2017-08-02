#!/usr/bin/env python

from __future__ import print_function

import argparse
import sys
import bed_functions

#FUNCTIONS#

def sort_end_pos(in_array, chroms):
    '''
    Takes in the data structure of SVs:
       [  [ ele0, ele1, ..., ele12, sv_max_dict, sv_95_dict ], ... ]
    Outputs the sorted index of the elements in an array eg. [5, 50, 32, 8, ..., n]
    Assumes b37 chromosmomes in the order of 1, 2, ..., 22, X, Y, MT, SVs involving other chromosomes are ignored
    TYPE:INTERCHROM is after all other types and is sorted by start_chrom then end_chrom then end_pos
    The ends of files are represented by #''# and are after all other SVs
    '''
    
    #Simple insertion sort#
    def SV_before(SV1, SV2):
        #Defines the sorting parameters. Roughly conforms to the default sorting of LUMPY output#
        #If SV2 is empty: return True
        if not SV2:
            return True
        elif not SV1:
            return False
        elif SV1[10] != 'TYPE:INTERCHROM' and SV2[10] == 'TYPE:INTERCHROM':
            return True
        elif SV1[10] == 'TYPE:INTERCHROM' and SV2[10] != 'TYPE:INTERCHROM':
            return False
        elif SV1[10] != 'TYPE:INTERCHROM': #Both not interchromosomal events
            if chroms[SV1[13]['end_chrom']] < chroms[SV2[13]['end_chrom']]:
                return True
            elif chroms[SV1[13]['end_chrom']] > chroms[SV2[13]['end_chrom']]:
                return False
            else: #On the same chromosome, if same end_pos same, arbitrarily return true
                if SV1[13]['end_pos'] <= SV2[13]['end_pos']:
                    return True
                else:
                    return False
        else: #Both are interchromosomal events
            if chroms[SV1[13]['start_chrom']] < chroms[SV2[13]['start_chrom']]:
                return True
            elif chroms[SV1[13]['start_chrom']] > chroms[SV2[13]['start_chrom']]:
                return False
            else: #Same starting chromosome
                if chroms[SV1[13]['end_chrom']] < chroms[SV2[13]['end_chrom']]:
                    return True
                elif chroms[SV1[13]['end_chrom']] > chroms[SV2[13]['end_chrom']]:
                    return False
                else: #Same ending chromosome; if same end_pos, arbitrarily return true
                    if SV1[13]['end_pos'] <= SV2[13]['end_pos']:
                        return True
                    else:
                        return False

    output = []
    for i in range(len(in_array)):
        if not in_array[i]:
            output.append(i)
        for j in range(len(output)):
            if SV_before(in_array[i], in_array[output[j]]):
                output.insert(j, i)
                break
        else:
            output.append(i)
    return output

def parse_line(line):
    '''
    Takes in a split input line and returns a SV data structure
    '''
    tmp1, tmp2 = line[13].split(';')
    sv_max = {}
    sv_95 = {}
    tmp1 = tmp1.split(':')
    tmp2 = tmp2.split(':')
    sv_max = { 'start_chrom': tmp1[1], 'start_pos': int(tmp1[2]), 'end_chrom': tmp2[0], 'end_pos': int(tmp2[1]) }
    tmp1, tmp2 = line[14].split(';')
    tmp1 = tmp1.split(':')
    tmp2 = tmp2.split(':')
    tmp3 = tmp1[2].split('-')
    tmp4 = tmp2[1].split('-')
    sv_95 = { 'start_chrom': tmp1[1], 'start_range': (int(tmp3[0]), int(tmp3[1])) , 'end_chrom': tmp2[0], 'end_range': (int(tmp4[0]), int(tmp4[1])) }
    evidence = []
    tmp1, tmp2  = line[11].split(':')
    if ';' in tmp2:
        evidence = ['split', 'discordant']
    else:
        if tmp2[0] == '1':
            evidence = ['discordant']
        else:
            evidence = ['split']    
    return line[:11] + [evidence] + line[12:13] + [sv_max, sv_95]
    
#PROGRAM#

parser = argparse.ArgumentParser(description='This program merges structural variants from multiple files of LUMPY output.')
parser.add_argument('-outfile', help="A file to write the merged SVs to")
parser.add_argument('infiles', nargs='+', help="Input files containing filtered output from LUMPY")
parser.add_argument('-bedfile', help="A bed file of exons to find overlapping SVs")

args = parser.parse_args()

FHs = []
current = []

#Assumes b37 chromosomes.
#Non cononical chromosomes are ignored
chroms = { '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9': 9, '10': 10, '11': 11, '12': 12,
           '13': 13, '14': 14, '15': 15, '16': 16, '17': 17, '18': 18, '19': 19, '20': 20, '21': 21, '22': 22,
           'X': 23, 'Y': 24, 'MT': 25 }

#Open all files simultaneously
for infile in args.infiles:
    FHs.append(open(infile, 'r'))

outfile = open(args.outfile, 'w')
bedfile = open(args.bedfile, 'r')

#Read first line of all files.
for fh in FHs:
    line = fh.readline().split()
    if not line:
        continue
    current.append(parse_line(line))

while True:
    #Look for non-cononical chromosomes. If found, read in a new SV#
    for i in range(len(current)):
        if not current[i]:
            continue
        x = True
        while x:
            if (not current[i][13]['start_chrom'] in chroms) or (not current[i][13]['end_chrom'] in chroms):
                line = FHs[i].readline().split()
                if not line:
                    current[i] = ''
                    x = False
                else:
                    current[i] = parse_line(line)
            else:
                x = False
    #Guarenteed to have SVs between cononical chromosomes or be at the end of a file#
    #Sort by end position#
    sorted_indecies = sort_end_pos(current, chroms)
    if not current[sorted_indecies[0]]: #If all files are done#
        sys.exit('Done')

#    print(str(sorted_indecies))
    #Start with the first index and try to merge it to the next index
    chr_start = ''
    chr_end = ''
    start_pos_front = ''
    start_pos_end = ''
    start_pos_max = ''
    end_pos_front = ''
    end_pos_end = ''
    end_pos_max = ''
    sv_type = ''
    start_uncertain = ''
    end_uncertain = ''
    evidence = [] #[split, discordant]
    ids = []
    ids_index = []
    for index in sorted_indecies:
#        print(str(current[index]))
        if not chr_start:
            chr_start = current[index][13]['start_chrom']
            chr_end = current[index][13]['end_chrom']
            start_pos_front = int(current[index][1])
            start_pos_end = int(current[index][2])
            start_pos_max = current[index][13]['start_pos']
            end_pos_front = int(current[index][4])
            end_pos_end = int(current[index][5])
            end_pos_max = current[index][13]['end_pos']
            sv_type = current[index][10][5:]
            start_uncertain = start_pos_end - start_pos_front
            end_uncertain = end_pos_end - end_pos_front
            ids.append(FHs[index].name)
            ids_index.append(index)
            evidence = current[index][11]
        else: #Not the first element
            if not current[index]:
                break
            elif (end_pos_end < int(current[index][4]) or 
                chr_end != current[index][13]['end_chrom'] or 
                chr_start != current[index][13]['start_chrom']): #Stop looking for merges if current SV is outside of the 95% range or is on a different chromosome#
                break
            else:
                if sv_type != current[index][10][5:]:
                    continue
                else: #Same type of event and same chromosome#
                    if (((start_pos_front <= int(current[index][1]) and start_pos_end >= int(current[index][1])) or
                         (start_pos_front <= int(current[index][2]) and start_pos_end >= int(current[index][2])) or
                         (start_pos_front <= int(current[index][1]) and start_pos_end >= int(current[index][2])) or
                         (start_pos_front >= int(current[index][1]) and start_pos_end <= int(current[index][2]))) and
                        ((end_pos_front <= int(current[index][4]) and end_pos_end >= int(current[index][4])) or
                         (end_pos_front <= int(current[index][5]) and end_pos_end >= int(current[index][5])) or
                         (end_pos_front <= int(current[index][4]) and end_pos_end >= int(current[index][5])) or
                         (end_pos_front >= int(current[index][4]) and end_pos_end <= int(current[index][5])))): #Merge SVs
                    #if (start_pos_front < current[index][13]['start_pos'] and start_pos_end > current[index][13]['start_pos'] and end_pos_front < current[index][13]['end_pos'] and end_pos_end > current[index][13]['end_pos']): #Merge SVs
                        if start_uncertain > (int(current[index][2]) - int(current[index][1])): #New has less uncertainty
                            start_pos_front = int(current[index][1])
                            start_pos_end = int(current[index][2])
                            start_pos_max = current[index][13]['start_pos']
                            start_uncertain = start_pos_end - start_pos_front
                        if end_uncertain > (int(current[index][5]) - int(current[index][4])):
                            end_pos_front = int(current[index][4])
                            end_pos_end = int(current[index][5])
                            end_pos_max = current[index][13]['end_pos']
                            end_uncertain = end_pos_end - end_pos_front
                        for evi in current[index][11]:
                            if not evi in evidence:
                                evidence.append(evi)
                        ids.append(FHs[index].name)
                        ids_index.append(index)

    #Relevant SVs are merged, print output#
    genes = []
    if chr_start == chr_end:
        genes = bed_functions.overlap(bedfile, chr_start, start_pos_max, end_pos_max)
        print('\t'.join([chr_start, str(start_pos_max), chr_end, str(end_pos_max), str(end_pos_max - start_pos_max), sv_type, ','.join(evidence), ','.join(ids), ','.join(genes)]), file=outfile)
        print('\t'.join([chr_start, str(start_pos_max), chr_end, str(end_pos_max), str(end_pos_max - start_pos_max), sv_type, ','.join(evidence), ','.join(ids), ','.join(genes)]))
    else:
        genes = bed_functions.overlap(bedfile, chr_start, start_pos_max - 1000, start_pos_max + 1000)
        genes = genes + bed_functions.overlap(bedfile, chr_end, end_pos_max - 1000, end_pos_max + 1000)
        print('\t'.join([chr_start, str(start_pos_max), chr_end, str(end_pos_max), '-', sv_type, ','.join(evidence), ','.join(ids), ','.join(genes)]), file=outfile)
        print('\t'.join([chr_start, str(start_pos_max), chr_end, str(end_pos_max), '-', sv_type, ','.join(evidence), ','.join(ids), ','.join(genes)]))

    #Read in new lines#
    for index in ids_index:
        line = FHs[index].readline().split()
        if not line:
            current[index] = ''
        else:
            current[index] = parse_line(line)
