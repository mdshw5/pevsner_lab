#!/usr/bin/env python3

import argparse
import sys
import math

def print_region(active, end, mle, sample, outfile):
    #active format:#
    #[sample number, rev_mle, seed_pos, start_pos, end_pos, forward_mle, current_mle]#
    out_arr = []
    out_arr.append('\t'.join(map(str, active[3])))
    out_arr.append('\t'.join(map(str, active[2])))
    out_arr.append('\t'.join(map(str, end)))
    out_arr.append(str(int(end[1]) - int(active[3][1])))
    out_arr.append(sample)
    out_arr.append(str(mle))
    print('\t'.join(out_arr), file=outfile)
                        
def print_header(outfile):
    print('start_chrom\tstart_pos\tseed_shrom\tseed_pos\tend_chrom\tend_pos\tregion_size\tsample_id\tmle', file=outfile)

parser = argparse.ArgumentParser(description='Detect mosaic variants from a VCF file')
parser.add_argument('invcf', metavar='.vcf', help='A VCF format file. Must have the "AD" tag in the genotype field')
parser.add_argument('outfile', metavar='.txt', help='An file to output potential mosaic site to')
parser.add_argument('--p_value', metavar='0.00-1.00', type=float, default='0.05', help='A p-value cutoff for outputting sites')
parser.add_argument('--seed', metavar='0.00-1.00', type=float, default='0.10', help='A probability cutoff for the initial detection of a mosiaic event.\
 Very low values will miss true positive events. High p-values will greatly increase runtime')
parser.add_argument('--depth', metavar='0+', type=int, default='10', help='A minimum depth for per sample filtering')
parser.add_argument('--max_prob', metavar='0.00-1.00', type=float, default='1e-3', help='Maximum probability of a single variant being mosaic or non-mosaic')
parser.add_argument('--filt_pass', action='store_true', help='Count only variants which are marked with "PASS" in the filter column')
parser.add_argument('--drop', metavar='0.00-1.00', type=float, default='0.1', help='During extension, if current score is (`--drop` * max_score) less than max_score, extension is stopped') 
args = parser.parse_args()

#Create bins for statistics#
histo = [0] * 1001
count = 0
total = 0
outfile = open(args.outfile, 'w')
print_header(outfile)
samples = []
with open(args.invcf, 'r') as f:
    ad = False
    for line in f:
        if line[0:15] == '##FORMAT=<ID=AD':
            ad = True
            continue
        elif line[0:2] == '##': #Skip header
            continue
        elif line[0:4] == '#CHR': #Populate sample list
            line = line.rstrip().split()
            for sample in line[9:]:
                samples.append(sample)
            continue
        line = line.rstrip().split()
        if args.filt_pass:
            if line[6] != 'PASS':
                continue
        if len(line[3]) > 1 or len(line[4]) > 1: #Skip indels and triallelic (+) sites
            continue
        if not ad:
            print('WARNING: Allele depth is not specified as a format tag.', file=sys.stderr)
        form = line[8].split(':')
        GT = 0
        AD = 0
        for i in range(len(form)): #Get genotype and allele depth from format field
            if form[i] == 'GT':
                GT = i
            elif form[i] == 'AD':
                AD = i
        for genotype in line[9:]: #For each sample
            genotype = genotype.split(':')
            if genotype[GT] == '0/1': #Ignore non-hets
                allele_counts = genotype[AD].split(',')
                if len(allele_counts) == 1: #Ignore sites without 2 alleles
                    continue
                ref_reads = allele_counts[0]
                alt_reads = allele_counts[1]
                if int(ref_reads) + int(alt_reads) < args.depth: #Sample depth filter
                    continue
                result = float(alt_reads) / (float(ref_reads) + float(alt_reads)) * 100
                #Add variant information#
                total += result
                count += 1
                result *= 10
                histo[int(round(result))] += 1

    #Print summaries#
#    print(str(histo), file=sys.stderr)
    average = total / float(count)
    print("Average is: " + str(average), file=sys.stderr)
    print("Total: " + str(total) + ' Count is: ' + str(count), file=sys.stderr)
    
    #Begin analysis#
    f.seek(0,0) #Seek to file start
    active = [] #Holds currently active mosaic events
    region = [] #Keeps track of chr and pos of previous variants
    prob_bin = False
    for line in f:
        if line[0] == '#':
            continue
        line = line.rstrip().split()
        if args.filt_pass:
            if line[6] != 'PASS':
                continue
        if not prob_bin: #Initalize prob_bin
            prob_bin = []
            for i in range(len(line[9:])):
                prob_bin.append([])
        if len(line[3]) > 1 or len(line[4]) > 1: #Skip indels or triallelic (+) sites
            continue
        form = line[8].split(':')
        GT = 0
        AD = 0
        if region:
            if region[-1][0] != line[0]:
                for j in range(len(active)):
                    if active[j][1] + active[j][5] < math.log(args.p_value):
                        print_region(active[j], active[j][4], active[j][1] + active[j][5], samples[active[j][0]], outfile)
                active = []
                region = []
                prob_bin = []
                for i in range(len(line[9:])):
                    prob_bin.append([])
        region.append((line[0], line[1]))
        for i in range(len(form)): #Get genotype and allele depth positions
            if form[i] == 'GT':
                GT = i
            elif form[i] == 'AD':
                AD = i
        for m in range(len(line[9:])): #m = sample index
            genotype = line[m + 9].split(':')
            if genotype[GT] == '0/1':
                allele_counts = genotype[AD].split(',')
                if len(allele_counts) == 1: #Skip sites w/o 2 alleles
                    print('INFO: Only one allele counted at: ' + str(line) + ' and index ' + str(m), file=sys.stderr)
                    prob_bin[m].append(('-', '-'))
                    continue
                ref_reads = allele_counts[0]
                alt_reads = allele_counts[1]
                if int(ref_reads) + int(alt_reads) < args.depth: #Sample depth filter
#                    print('INFO: Depth cutoff not met at: ' + str(line) + ' and index ' + str(m), file=sys.stderr)
                    prob_bin[m].append(('-', '-'))
                    continue
                result = float(alt_reads) / (float(ref_reads) + float(alt_reads)) * 100
                #p_out = probablity variant does not fall in normal distribution#
                #p_in = probability variant falls in normal distribution#
                if result > average:
                    p_out = (sum(histo[int(round(result * 10)):]) + sum(histo[:int(round((average + (average - result)) * 10))])) / float(count)
                    p_in = 1 - p_out
                else:
                    p_out = (sum(histo[:int(round(result * 10))]) + sum(histo[int(round((average + (average - result)) * 10)):])) / float(count)
                    p_in = 1 - p_out
                #Scale by max/min probabilities#
                if p_out < args.max_prob:
                    p_out = args.max_prob
                    p_in = 1 - args.max_prob
                elif p_in < args.max_prob:
                    p_in = args.max_prob
                    p_out = 1 - args.max_prob
                #Convert to natural log probabilities#
                log_p_out = math.log(p_out)
                log_p_in = math.log(p_in)
                if m in [j[0] for j in active]: #If the sample is marked as active
                    #Forward Extend#
                    j = [p[0] for p in active].index(m)
                    #j = Index of the sample in the `active` list#
                    if active[j][2][0] != region[-1][0]: #Do not extend across chromosomes
                        if active[j][1] + active[j][5] < math.log(args.p_value):
                            print_region(active[j], active[j][4], active[j][1] + active[j][5], samples[m], outfile)
                        active.pop(j)
                    elif active[j][2] in region: #The seed can be found
#                        k = region.index(active[j][2])
                        active[j][6] += log_p_out - log_p_in
                        if active[j][6] > 0 or active[j][6] > active[j][5] - (active[j][5] * args.drop): #Stop forward extension
                            if active[j][1] + active[j][5] < math.log(args.p_value):
                                print_region(active[j], active[j][4], active[j][1] + active[j][5], samples[m], outfile)
                            active.pop(j)
                        elif active[j][6] < active[j][5]: #Current is better than best
                            active[j][5] = active[j][6]
                            active[j][4] = region[-1]                        
                    else:
                        print("WARNING: Region start not found with: " + str(active[j]), file=sys.stderr)
                        if active[j][1] + active[j][5] < math.log(args.p_value):
                            print_region(active[j], active[j][4], active[j][1] + active[j][5], samples[m], outfile)
                        active.pop(j)
                elif p_out < args.seed:
                    #Reverse Extend
                    mle = log_p_out - log_p_in
                    current = log_p_out - log_p_in
                    start = region[-1]
                    for j in range(len(prob_bin[i]) -2, -1, -1):
                        #j=variant index in region + 1 or index in prob_bin
                        if region[j - 1][0] != region[-1][0]:
                            break
                        if prob_bin[m][j][0] == '-':
                            continue
                        current += prob_bin[m][j][0] - prob_bin[m][j][1]
                        if current < mle:
                            mle = current
                            start = region[j - 1]
                        if current > 0 or current > mle - (mle * args.drop):
                            break
                    #Correct mle for addition of seed to start and end extension#
                    mle = mle - (log_p_out - log_p_in)
                    active.append([m, mle, region[-1], start, region[-1], log_p_out - log_p_in, log_p_out - log_p_in]) #i (sample number), rev_mle,                                                                                                          #seed pos, start pos, #end pos
                                                                                                           #fwd_mle, fwd_curr_mle 
                prob_bin[m].append((log_p_out, log_p_in))
            else:
                prob_bin[m].append(('-', '-'))

                
