#!/usr/bin/env bash
#
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 7-00:05
#SBATCH --mem 4096

# $1 = VCF
# $2 = REF

java -jar $BCBIOVAR variant-prep $1 $2