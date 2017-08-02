#!/usr/bin/env bash

gzip -dc $1 | sambamba view -S -f bam -l 0 /dev/stdin | sambamba sort --tmpdir /mnt/data/don/Tmp -m 60GB -t 16 -o $2 /dev/stdin