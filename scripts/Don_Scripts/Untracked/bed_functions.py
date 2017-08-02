#!/usr/bin/env python

from __future__ import print_function

def overlap(infile, chrom, start, stop):
    '''
    Input is an open filehandle of a sorted bed file, chromosome, start pos and stop pos.
    Output is a list of elements from the bed file overlapping those coordinates
    '''
    genes = []
    reset = False
    before_start = False
    same_chrom = False
    while True:
        line = infile.readline()
        if (not line):
            if reset:
                return genes
            elif not genes:
                reset = True
                before_start = True
                infile.seek(0, 0)
                line = infile.readline()
            else:
                return genes
        line = line.split()
#        print(str(line))
        if line[0] == chrom:
            same_chrom = True
            if before_start:
                if int(line[1]) > int(stop):
                    return genes
                else:
                    if int(line[1]) >= int(start):
                        if not line[3] in genes:
                            genes.append(line[3])
                    else:
                        if int(line[2]) >= int(start):
                            if not line[3] in genes:
                                genes.append(line[3])
            else:
                if int(line[1]) < int(start):
                    before_start = True
                    if int(line[2]) > int(stop):
                        if not line[3] in genes:
                            genes.append(line[3])
                else:
                    pnt_pos = infile.tell()
                    if pnt_pos <= 2000:
                        before_start = True
                        reset = True
                        infile.seek(0, 0)
                    else:
                        infile.seek(-2000, 1)
                        infile.readline()
        else:
            if same_chrom:
                before_start = True
                continue
            elif reset:
                continue
            else:
                reset = True
                before_start = True
                infile.seek(0, 0)

def nearest_upstream(infile, chrom, pos):
    reset = False
    same_chrom = False
    before_start = False
    last = ''
    while True:
        line = infile.readline()
        if not line:
            if not last:
                reset = True
                before_start = True
                infile.seek(0, 0)
            else:
                return last
        line = line.split()
        if line[0] == chrom:
            same_chrom = True
            if before_start:
                if int(line[2]) < int(pos):
                    last = line[3]
                else:
                    if int(line[1]) > int(pos):
                        return last
            else:
                if int(line[2]) < int(pos):
                    before_start = True
                    last = line[3]
                else:
                    pnt_pos = infile.tell()
                    if pnt_pos <= 2000:
                        before_start = True
                        reset = True
                        infile.seek(0, 0)
                    else:
                        infile.seek(-2000, 1)
                        infile.readline()
        else:
            if same_chrom:
                before_start = True
                continue
            elif reset:
                continue
            else:
                reset = True
                before_start = True
                infile.seek(0, 0)

def nearest_downstream(infile, chrom, pos):
    reset = False
    same_chrom = False
    before_start = False
    while True:
        line = infile.readline()
        if not line:
            if not reset:
                reset = True
                before_start = True
                infile.seek(0, 0)
            else:
                return ''
        else:
            return ''
        if line[0] == chrom:
            same_chrom = True
            if before_start:
                if int(pos) < int(line[1]):
                    return line[3]
            else:
                if int(pos) > int(line[2]):
                    before_start = True
                else:
                    if int(pos) < int(line[1]):
                        pnt_pos = infile.tell()
                        if pnt_pos <= 2000:
                            before_start = True
                            reset = True
                            infile.seek(0, 0)
                        else:
                            infile.seek(-2000, 1)
                            infile.readline()
        else:
            if same_chrom:
                before_start = True
                continue
            elif reset:
                continue
            else:
                reset = True
                before_start = True
                infile.seek(0, 0)

if __name__ == "__main__":
    import sys
    f = open(sys.argv[1], 'r')
    print(str(overlap(f, sys.argv[2], sys.argv[3], sys.argv[4])))
