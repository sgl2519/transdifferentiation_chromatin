#!/bin/bash
# get list of non-redundant TSS for each gene contained in the gtf file 
# remember to filter out non protein-coding genes from the gtf in case you don't want them

gtf_file="$1"

classified_exons=$(gawk -f ./classify_exons.sh $gtf_file)
for strand in - +; do
	if [ "$strand" == "-" ]; then
		col=3
	else
		col=2
	fi
	grep -F "$strand" <(grep "unique\|first" <(echo "$classified_exons") | awk 'BEGIN{FS=OFS="\t"}{print $1,($4-1),$5,$11,"0","-",$10}') | sort -u -k7,7 -k"$col","$col"n;
done


