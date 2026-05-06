#!/bin/bash

# Promoter marks

for i in H3K4me1 H3K4me2 H3K4me3 H3K27ac H3K9ac H3K9me3 H3K27me3

do

cd ../analysis/all.marks/"$i"/

#*****
# 1. prepare folder 
# to store results of 
# intersecting peaks
mkdir peaks 
#*****

#*****
# 2. Intersect Zerone peaks 
# against the genomic interval
# [1st region of interest (-2 Kb of 1st TSS), last TTS] 
# all time-points

# 2.1. define myRange
# (2000 for all marks except K36me3, K20me1 --> 0)
myRange=2000 

# 2.2. run intersection
cat "$i".peaks.txt | while read file; do bedtools intersect -a <(awk -v myRange=$myRange 'BEGIN{FS=OFS="\t"}{if ($6=="+") {$2=$2-myRange} else {$3=$3+myRange}; print $0}' ../../../../bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed) -b $file > peaks/"$(echo $file | cut -d/ -f 10)".intersectBed.bed; done
#*****

cd QN.merged

#*****
# 3. add link to gene expression metadata within the folder
ln -s ../expression/QN.merged/metadata.class2.tsv metadata.tsv

#*****
# 4. compute peak length per gene in the different time points + gene length

# 4.1. peak length at 0h
#conda deactivate
../../../../bin/join.py -b <(tail -n+2 metadata.tsv | cut -f1) -a <(cat ../peaks/H000.intersectBed.bed | awk '{print ($3-$2)"\t"$4}'| cut -d "." -f1 | awk '{print $2"\t"$1}' | awk '{a[$1]+=$2}END{for(k in a)print k"\t"a[k]}') -u -p 0 | sed '1iH000' > all.tp.peak.length.QN.merged.tsv 

# 4.2. peak length at other 11 time-points
for time in H003 H006 H009 H012 H018 H024 H036 H048 H072 H120 H168; do ../../../../bin/join.py -b all.tp.peak.length.QN.merged.tsv -a <(cat ../peaks/"$time".intersectBed.bed | awk '{print ($3-$2)"\t"$4}'| cut -d "." -f1 | awk '{print $2"\t"$1}' | awk '{a[$1]+=$2}END{for(k in a)print k"\t"a[k]}') -u -p 0 --b_header | sed "1s/V1/$time/" > tmp; mv tmp all.tp.peak.length.QN.merged.tsv; done 

# 4.3. add range of extension to gene body
# (0 --> K36me3, K20me1; 2000 --> all other marks)
myRange=2000 

# 4.4. compute gene length
../../../../bin/join.py -b all.tp.peak.length.QN.merged.tsv -a <(awk -v myRange=$myRange 'BEGIN{FS=OFS="\t"}{if ($6=="+") {$2=$2-myRange} else {$3=$3+myRange}; print $0}' ../../../../bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed | awk '{print $3-$2"\t"$4}' | cut -d "." -f1 | awk '{print $2"\t"$1}') --b_header | sed '1s/V1/gene_length/' > tmp; mv tmp all.tp.peak.length.QN.merged.tsv 
#*****

#*****
# 6. select genes with a peak in the region of interest in at least one time point

# 6.1. define regions of interest
# gene body --> K36me3, K20me1
# [-2 Kb, +2 Kb] --> all other marks
myBed="../../../../bed.files/gencode.v24.protein.coding.non.overlapping.genes.5Kb.upstream.downstream.TSS.bed" 

# 6.2. compute intersection
cat ../"$i".peaks.txt | while read file; do bedtools intersect -a $myBed -b $file; done | cut -f7 | sort -u | cut -d "." -f1 > all.genes.intersecting.peaks.tsv 
grep -Fx -f all.genes.intersecting.peaks.tsv <(tail -n+2 expression.matrix.tsv | cut -f1) > genes.intersecting.peaks.tsv
#*****


done


# Gene body marks

for i in H3K36me3 H4K20me1

do

cd ../../../../analysis/all.marks/"$i"/

#*****
# 1. prepare folder 
# to store results of 
# intersecting peaks
mkdir peaks 
#*****

#*****
# 2. Intersect Zerone peaks 
# against the genomic interval
# [1st region of interest (-2 Kb of 1st TSS), last TTS] 
# all time-points

# 2.1. define myRange
# (2000 for all marks except K36me3, K20me1 --> 0)
myRange=0 

# 2.2. run intersection
cat "$i".peaks.txt | while read file; do bedtools intersect -a <(awk -v myRange=$myRange 'BEGIN{FS=OFS="\t"}{if ($6=="+") {$2=$2-myRange} else {$3=$3+myRange}; print $0}' ../../../../bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed) -b $file > peaks/"$(echo $file | cut -d/ -f 10)".intersectBed.bed; done
#*****

cd QN.merged

#*****
# 3. add link to gene expression metadata within the folder
ln -s /no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/expression/QN.merged/metadata.class2.tsv metadata.tsv

#*****
# 4. compute peak length per gene in the different time points + gene length

# 4.1. peak length at 0h
#conda deactivate
../../../../bin/join.py -b <(tail -n+2 metadata.tsv | cut -f1) -a <(cat ../peaks/H000.intersectBed.bed | awk '{print ($3-$2)"\t"$4}'| cut -d "." -f1 | awk '{print $2"\t"$1}' | awk '{a[$1]+=$2}END{for(k in a)print k"\t"a[k]}') -u -p 0 | sed '1iH000' > all.tp.peak.length.QN.merged.tsv 

# 4.2. peak length at other 11 time-points
for time in H003 H006 H009 H012 H018 H024 H036 H048 H072 H120 H168; do /no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/bin/join.py -b all.tp.peak.length.QN.merged.tsv -a <(cat ../peaks/"$time".intersectBed.bed | awk '{print ($3-$2)"\t"$4}'| cut -d "." -f1 | awk '{print $2"\t"$1}' | awk '{a[$1]+=$2}END{for(k in a)print k"\t"a[k]}') -u -p 0 --b_header | sed "1s/V1/$time/" > tmp; mv tmp all.tp.peak.length.QN.merged.tsv; done 

# 4.3. add range of extension to gene body
# (0 --> K36me3, K20me1; 2000 --> all other marks)
myRange=0 

# 4.4. compute gene length
../../../../bin/join.py -b all.tp.peak.length.QN.merged.tsv -a <(awk -v myRange=$myRange 'BEGIN{FS=OFS="\t"}{if ($6=="+") {$2=$2-myRange} else {$3=$3+myRange}; print $0}' ../../../../bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed | awk '{print $3-$2"\t"$4}' | cut -d "." -f1 | awk '{print $2"\t"$1}') --b_header | sed '1s/V1/gene_length/' > tmp; mv tmp all.tp.peak.length.QN.merged.tsv 
#*****

#*****
# 6. select genes with a peak in the region of interest in at least one time point

# 6.1. define regions of interest
# gene body --> K36me3, K20me1
# [-2 Kb, +2 Kb] --> all other marks
myBed="../../../../bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed" 

# 6.2. compute intersection
cat ../"$i".peaks.txt | while read file; do bedtools intersect -a $myBed -b $file; done | cut -f7 | sort -u | cut -d "." -f1 > all.genes.intersecting.peaks.tsv 
grep -Fx -f all.genes.intersecting.peaks.tsv <(tail -n+2 expression.matrix.tsv | cut -f1) > genes.intersecting.peaks.tsv
#*****


done
