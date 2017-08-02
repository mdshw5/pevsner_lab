#!/usr/bin/env python3

import math
import shlex
import subprocess
import os.path
import sys
import os
import time

def genotype_bams(sample_dict, out_dir, reference, gatk, tmp_dir, memory='1g', \
                      dispatch='srun', max_jobs=50, intervals=50, trailing='_sorted.bam', \
                      args='--emitRefConfidence GVCF --variant_index_type LINEAR ' +\
                      '--variant_index_parameter 128000 -A BaseQualityRankSumTest -A RMSMappingQuality ' +\
                      '-A TandemRepeatAnnotator -A QualByDepth -A MappingQualityRankSumTest ' +\
                      '-A ReadPosRankSumTest -A FisherStrand -A Coverage -A MappingQualityZero ' +\
                      '-A SpanningDeletions -A StrandOddsRatio'):
    
    if not os.path.isfile(reference):
        raise ValueError("Reference file not found")
    if not os.path.isfile(reference + '.fai'):
        raise ValueError("Reference index file not found")
    for sample, bam_list in sample_dict.items():
        for bam in bam_list:
            if not os.path.isfile(bam):
                raise ValueError("BAM file not found: " + bam)
    
    total_size = 0
    chroms = []
    bp = []

    # Calculate the step size #
    with open(reference + '.fai') as f:
        for line in f:
            line = line.rstrip().split()
            total_size += int(line[1])
            chroms.append(line[0])
            bp.append(int(line[1]))

    step_size = math.ceil(total_size / intervals)
    
    # Get the intervals #
    intervals = []
    last_chrom = ''
    last_pos = 1
    passed = 0
    cur_interval = []
    
    for i in range(len(chroms)):
        if chroms[i] == 'hs37d5': # Skip bait chromosome
            continue
        while ( last_pos + step_size - passed < bp[i]):
            cur_interval.append('-L ' + chroms[i] + ':' + str(last_pos) + '-' + str(last_pos + step_size - passed))
            intervals.append(' '.join(cur_interval))
            cur_interval = []
            last_pos += step_size - passed
            passed = 0
        else:
            cur_interval.append('-L ' + chroms[i] + ':' + str(last_pos) + '-' + str(bp[i]))
            passed += bp[i] - last_pos
            last_pos = 1
    if cur_interval:
        intervals.append(' '.join(cur_interval))

    # Call variants over the intervals #
    running = 0
    jobs_queue = []
    bam_idx = 0
    variant_files = {}

    for sample, bam_list in sample_dict.items():
        for idx, interval in enumerate(intervals):
            infile = ' -I ' + ' -I '.join(bam_list)
            
            command = dispatch + ' -c 2 java -Xmx' + memory + ' -Djava.io.tmpdir=' + \
                tmp_dir + ' -jar ' + gatk + ' -T HaplotypeCaller -R ' + reference + infile + ' -o ' + tmp_dir + \
                '/' + sample + str(idx) + '.vcf ' + interval + ' ' + args

            # Check number of currently running jobs #
            while len(jobs_queue) >= max_jobs:
                time.sleep(20)
                to_pop = []
                for i in range(len(jobs_queue)):
                    if jobs_queue[i].poll() != None:
                        to_pop.append(i)
                for i in reversed(to_pop):
                    jobs_queue.pop(i)

            print("Running: " + command, file=sys.stderr)
            if not sample in variant_files:
                variant_files[sample] = []
            variant_files[sample].append(tmp_dir + '/' + sample + str(idx) + '.vcf')
            shell = shlex.split(command)
            jobs_queue.append(subprocess.Popen(shell))
            time.sleep(20)
    for job in jobs_queue:
        job.wait()
            
    # Combine the variants from a single individual #
    max_jobs = math.ceil(max_jobs / 10) # I/O bound
    jobs_queue = []
    gvcfs = []
    
    for name, vcf_list in variant_files.items():
        variants = ' -V ' + ' -V '.join(vcf_list)

        command = dispatch + ' -c 1 java -Xmx' + memory + ' -cp ' + gatk + \
            ' org.broadinstitute.gatk.tools.CatVariants -R ' + reference + variants + ' -out ' + \
            out_dir + '/' + name + '.g.vcf -assumeSorted'
    
        # Check the number of running jobs #
        while len(jobs_queue) >= max_jobs:
            time.sleep(20)
            to_pop = []
            for i in range(len(jobs_queue)):
                if jobs_queue[i].poll() != None:
                    to_pop.append(i)
            for i in reversed(to_pop):
                jobs_queue.pop(i)

        print("Running: " + command, file=sys.stderr)
        shell = shlex.split(command)
        jobs_queue.append(subprocess.Popen(shell))
        gvcfs.append(out_dir + '/' + name + '.g.vcf')
    for job in jobs_queue:
        job.wait()
        
    for name, vcf_list in variant_files.items():
        for vcf in vcf_list:
            if os.path.isfile(vcf):
                os.remove(vcf)
            if os.path.isfile(vcf + '.idx'):
                os.remove(vcf + '.idx')

    return gvcfs
