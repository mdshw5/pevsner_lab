#!/usr/bin/env bash

for infile in "$@"
do
    samtools index $infile
done