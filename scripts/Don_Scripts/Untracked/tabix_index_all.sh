#!/usr/bin/env bash

for infile in "$@"
do
    tabix -p vcf $infile
done