#!/usr/bin/env python3

import argparse
import gzip
import os
import re
import sys

def main(target, library = None, platform = "ILLUMINA", outdir = "data/processed_fastq"):

    if not (target.endswith('R1.fastq.gz') or target.endswith('R1.fastq')):
        sys.exit('File name must end in R1.fastq[.gz].')
    
    r1in = os.path.abspath(target)
    r2in = re.sub('_R1\.fastq', '_R2.fastq', r1in)
    print("Processing" + r1in + " and " + r2in)
    
    if not os.path.isfile(r2in):
        err = "R2 of " + target + " doesn't exist."
        sys.exit(err)

    sample_name = re.sub('.*/(.*)_R[1-2]{1}.[^/]*', '\\1', os.path.abspath(target))  #Adjusted for KK_700_series

    if target.endswith('gz'):
        f = gzip.open(target, 'rb')
    else:
        f = open(target, 'r')
    header = f.readline()

    if target.endswith('gz'):
        seq_ids = header.decode('UTF-8').rstrip().replace('@','').split(':')
    else:
        seq_ids = header.rstrip().replace('@','').split(':')

    seq_dict = dict()
    seq_dict['machine'] = seq_ids[0]
    seq_dict['readid'] = seq_ids[1]
    seq_dict['flowcell'] = seq_ids[2]
    seq_dict['lane'] = seq_ids[3]
    seq_dict['index'] = seq_ids[9]

    if library == None:
        lb = '-'.join([sample_name, seq_dict['index']])
    else:
        lb = library

    # symbolic link name: <sample_name>_<library>_<machine>_<flowcell_id>_<sample_barcode>_<lane>_R[1/2].fastq[.gz]
    sym_name = '_'.join([sample_name, lb, seq_dict['machine'], seq_dict['flowcell'], seq_dict['index'], seq_dict['lane']])

    if target.endswith('gz'):
        r1out = outdir + '/' + sym_name + '_R1.fastq.gz'
        r2out = outdir + '/' + sym_name + '_R2.fastq.gz'
    else:
        r1out = outdir + '/' + sym_name + '_R1.fastq'
        r2out = outdir + '/' + sym_name + '_R2.fastq'

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    os.symlink(r1in, r1out)
    os.symlink(r2in, r2out)

if __name__ == "__main__":

    def parse_arguments():
        parser = argparse.ArgumentParser(prog='PROG',
                                         description='Parse readgroups from FASTQ file')
        parser.add_argument('target', type=str, help='Target fastq[.gz] file to parse read groups')
        parser.add_argument('-lb', '--library', type=str, help='Library of the sample')
        parser.add_argument('-pl', '--platform', type=str, default='ILLUMINA', help='Sequencing platform, defaults to Illumina')
        parser.add_argument('-o', '--outdir', type=str, default='data/processed_fastq', help='Output directory for coerced symbolic links')
        return parser.parse_args()

    args = parse_arguments()

    main(target = args.target, library = args.library, platform = args.platform, outdir = args.outdir)
