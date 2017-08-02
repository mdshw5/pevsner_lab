#!/usr/bin/env python

from __future__ import print_function

def SV_before(SV1, SV2):
        #Defines the sorting parameters. Roughly conforms to the default sorting of LUMPY output#
        #If SV2 is empty: return True
    chroms = { '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6, '7': 7, '8': 8, '9':9, '10': 10, '11': 11, '12': 12,
           '13': 13, '14': 14, '15': 15, '16': 16, '17': 17, '18': 18, '19': 19, '20': 20, '21': 21, '22': 22,
           'X': 23, 'Y': 24, 'MT': 25 }
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

def parse_lumpy_line(line):
    '''
    Input is a line of LUMPY output split into an array
    Returns a SV data structure
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
