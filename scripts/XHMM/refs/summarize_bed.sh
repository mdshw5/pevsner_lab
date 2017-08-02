for b in *.bed
do
    cov=$(cat $b | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}')
    echo -e "$b\t$cov" >> bed_summaries.txt
done

