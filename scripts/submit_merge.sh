#!/usr/bin/env bash

a=()
d=$(sed 's/-//g' <<< $(date +%F))

for fq in data/fastq/*.gz
do
	foo=${fq%_R*}
	sample=${foo##*/}
	a+=($sample)
done

# Get unique samples up to Lane
u=(`echo "${a[@]}" | tr ' ' '\n' | sort -u | tr '\n' ' '`)

# Submit
for i in "${u[@]}"
do
    echo "Merging sample $i"
    sbatch scripts/do_merge.sh $i --error "logs/merge_fastq/${d}_${i}_merge_fastq.err" --output "logs/merge_fastq/${d}_${i}_merge_fastq.log"
done

# Document
dt=$(sed 's/-//g'<<< $(date "+%F [%T]"))
echo "${dt} : Merge sample_lane sets [${0}]" >> CHANGELOG
