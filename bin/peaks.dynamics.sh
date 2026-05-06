#*****
# get peaks dynamics

#!/bin/bash

# Promoter regions

for i in H3K4me1 H3K4me2 H3K4me3 H3K27ac H3K9ac H3K9me3 H3K27me3

do

cd ../analysis/all.marks/"$i"/QN.merged/

mkdir ant.del.analysis

cd ant.del.analysis

myBed="../../../../../bed.files/gencode.v24.protein.coding.non.overlapping.genes.5Kb.upstream.downstream.TSS.bed" &&
cat ../../"$i".peaks.txt | while read file; do bedtools intersect -a $myBed -b $file | awk -v file=$(echo $file | cut -d/ -f 10) '{print $0"\t"file}' | cut -f7-8 | sort -u; done | /users/rg/bborsari/bin/make.column.list.py > "$i".peaks.dynamics.tsv

done


# Gene body regions

for i in H3K36me3 H4K20me1

do

cd ../../../../../analysis/all.marks/"$i"/QN.merged/

mkdir ant.del.analysis

cd ant.del.analysis

myBed="../../../../../bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed" &&
cat ../../"$i".peaks.txt | while read file; do bedtools intersect -a $myBed -b $file | awk -v file=$(echo $file | cut -d/ -f 10) '{print $0"\t"file}' | cut -f7-8 | sort -u; done | /users/rg/bborsari/bin/make.column.list.py > "$i".peaks.dynamics.tsv

done
