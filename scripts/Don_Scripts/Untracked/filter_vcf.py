#!/usr/bin/env python

from __future__ import print_function

import argparse
import sys

parser = argparse.ArgumentParser(description='Filter a vcf based on parameters in a vcf. Filters will exclude (filter) values equal to the filter')
parser.add_argument('-BQRS', type=float, help='Symmetric filter of base quality score rank sum test')
parser.add_argument('-CLPRS', type=float, help='Symmetric filter of clipping rank sum test')
parser.add_argument('-depth', type=int, help='Filter by a depth requiremnet')
parser.add_argument('-FSBias', type=float, help='Symmetric filter of the result of a fisher test for strand bias')
parser.add_argument('-mapq', type=float, help='Filter by average mapq')
parser.add_argument('-MQRS', type=float, help='Symmetric filer of mapping quality rank sum test')
parser.add_argument('-qual', type=float, help='Filter by QUAL field')
parser.add_argument('-ReadPosRS', type=float, help='Symmetric filter of read position rank sum test')
parser.add_argument('-Inbreeding', type=float, help='Symmetric filter of inbreeding coefficient')
parser.add_argument('-MQ0', type=int, help='Maximum number of mapping quality zero reads')
parser.add_argument('-QD', type=float, help='Minimum variant quality by depth')
parser.add_argument('vcf', help='Input vcf format file')
parser.add_argument('filtered', help='Output filtered vcf file')

args = parser.parse_args()

outfile = open(args.filtered, 'w')
with open(args.vcf, 'r') as f:
    for line in f:
        line = line.rstrip()
        if line[0] == '#':
            print(line, file=outfile)
        else:
            line = line.split()
            if args.qual:
                if float(line[5]) <= args.qual:
                    line[6] = 'QUAL'
                    print('\t'.join(line), file=outfile)
                    continue
            info_tags = line[7].split(';')
            info_dict = {}
            for tag in info_tags:
                tmp = tag.split('=')
                if len(tmp) == 1:
                    info_dict[tmp[0]] = ''
                else:
                    info_dict[tmp[0]] = tmp[1]
            if args.BQRS:
                if 'BaseQRankSum' in info_dict:
                    if (float(info_dict['BaseQRankSum']) <= (-args.BQRS) or float(info_dict['BaseQRankSum']) >= args.BQRS):
                        line[6] = 'BaseQRankSum'
                        print('\t'.join(line), file=outfile)
                        continue
                else:
                    print("No BQRS at chromosome " + line[0] + " position " + line[1], file=sys.stderr)
            if args.CLPRS:
                if 'ClippingRankSum' in info_dict:
                    if (float(info_dict['ClippingRankSum']) <= (-args.CLPRS) or float(info_dict['ClippingRankSum']) >= args.CLPRS):
                        line[6] = 'ClippingRankSum'
                        print('\t'.join(line), file=outfile)
                        continue
                else:
                    print("No ClippingRankSum at chromosome " + line[0] + " position " + line[1], file=sys.stderr)
            if args.depth:
                if 'DP' in info_dict:
                    if (int(info_dict['DP']) <= args.depth):
                        line[6] = 'Depth'
                        print('\t'.join(line), file=outfile)
                        continue
                else:
                    print("No DP at chromosome " + line[0] + " position " + line[1], file=sys.stderr)
            if args.FSBias:
                if 'FS' in info_dict:
                    if (float(info_dict['FS']) <= (-args.FSBias) or float(info_dict['FS']) >= args.FSBias):
                        line[6] = 'FS'
                        print('\t'.join(line), file=outfile)
                        continue
                else:
                    print("No FS at chromosome " + line[0] + " position " + line[1], file=sys.stderr)
            if args.mapq:
                if 'MQ' in info_dict:
                    if (float(info_dict['MQ']) <= args.mapq):
                        line[6] = 'MQ'
                        print('\t'.join(line), file=outfile)
                        continue
                else:
                    print("No MQ at chromosome " + line[0] + " position " + line[1], file=sys.stderr)
            if args.MQRS:
                if 'MQRankSum' in info_dict:
                    if (float(info_dict['MQRankSum']) <= (-args.MQRS) or float(info_dict['MQRankSum']) >= args.MQRS):
                        line[6] = 'MQRankSum'
                        print('\t'.join(line), file=outfile)
                        continue
                else:
                    print("No MQRankSum at chromosome " + line[0] + " position " + line[1], file=sys.stderr)
            if args.ReadPosRS:
                if 'ReadPosRankSum' in info_dict:
                    if (float(info_dict['ReadPosRankSum']) <= (-args.ReadPosRS) or float(info_dict['ReadPosRankSum']) >= args.ReadPosRS):
                        line[6] = 'ReadPosRankSum'
                        print('\t'.join(line), file=outfile)
                        continue
                else:
                    print("No ReadPosRankSum at chromosome " + line[0] + " position " + line[1], file=sys.stderr)
            if args.Inbreeding:
                if 'InbreedingCoeff' in info_dict:
                    if (float(info_dict['InbreedingCoeff']) <= (-args.Inbreeding) or float(info_dict['InbreedingCoeff']) >= args.Inbreeding):
                        line[6] = 'InbreedingCoeff'
                        print('\t'.join(line), file=outfile)
                        continue
                else:
                    print("No InbreedingCoeff at chromosome " + line[0] + " position " + line[1], file=sys.stderr)
            if args.MQ0:
                if 'MQ0' in info_dict:
                    if (int(info_dict['MQ0']) <= args.MQ0):
                        line[6] = 'MQ0'
                        print('\t'.join(line), file=outfile)
                        continue
                else:
                    print("No MQ0 at chromosome " + line[0] + " position " + line[1], file=sys.stderr)
            if args.QD:
                if 'QD' in info_dict:
                    if (float(info_dict['QD']) <= args.QD):
                        line[6] = 'QD'
                        print('\t'.join(line), file=outfile)
                        continue
                else:
                    print("No QD at chromosome " + line[0] + " position " + line[1], file=sys.stderr)
            line[6] = 'PASS'
            print('\t'.join(line), file=outfile)
