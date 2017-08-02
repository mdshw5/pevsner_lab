#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='Take VAAST .simple output and transform it to excel input format')
parser.add_argument('infile', help='A VAAST .simple file for input.')
args = parser.parse_args()

with open(args.infile) as f:
    for line in f:
        if line.startswith('RANK'):
            print('RANK\tGENE\tP-value\tScore\tVar1\tVar1_score\tVar1_Nucleotide_change\tVar1_AA_change\tVar1_background_alleles\tVar1_case_alleles\tVar2\tVar2_score\tVar2_Nucleotide_change\tVar2_AA_change\tVar2_background_alleles\tVar2_case_alleles\t...')
            continue
        else:
            line = line.rstrip().split()
            if float(line[4]) == 0:
                continue
            outline = line[0:3]
            outline.append(line[4])
            for variant in line[5:]:
                variant = variant.split(';')
                outline.append(variant[0])
                outline.append(variant[1])
                outline.append(variant[2])
                outline.append(variant[3])
                ref_allele, alt_allele = variant[4].split(',')
                outline.append(ref_allele)
                outline.append(alt_allele)
            print('\t'.join(outline))
