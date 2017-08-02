#!/usr/bin/env bash

for infile in "$@"
do
    /mnt/data/don/Scripts/Tracked/sort_lumpy_sv.py --infile $infile > ${infile%.lumpy_out.txt}.lumpy_sorted.txt
done