#!/bin/bash

for i in H3K36me3 H4K20me1

do

## Start process

#******
# 1. prepare folders for logs and errors
mkdir ../analysis/all.marks/"$i"/
cd ../all.marks/"$i"/
mkdir logs errors
#******

#******
# 3. define bedfile for this mark
myBed="../../bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed" 
#******

#******
# 4. list pileup bw files
grep "${i}X1" ../../../references/pipeline.db | sort -k1,1 | cut -f2 > "$i".R1.pvalue.bw.txt 
grep "${i}X2" ../../../references/pipeline.db | sort -k1,1 | cut -f2 > "$i".R2.pvalue.bw.txt 
#******

#******
# 5. list Zerone peak calling files
ls folder/to/bed/files | sort -k1 > "$i".peaks.txt 
#******

#******
# 6. prepare folder for bwtool analysis
mkdir bwtool 
#******

#******
# 7. get mark matrix

# 7.1. rep. 1 
../../../bin/get.matrix.chipseq.sh --bedfile $myBed --bw "$i".R1.pvalue.bw.txt --target "$i" --outFolder bwtool --signal mean --peaks "$i".peaks.txt --keep yes --outFile "$i".R1.matrix.tsv 

# 7.2. clean-up
mkdir coordinates.bed 
mv *.intersection.bed coordinates.bed/ 
rm "$i".path.* 

# 7.3. rep. 2
../../../bin/get.matrix.chipseq.sh --bedfile $myBed --bw "$i".R2.pvalue.bw.txt --target "$i" --outFolder bwtool --signal mean --peaks "$i".peaks.txt --keep yes --outFile "$i".R2.matrix.tsv 

# 7.4. clean-up
rm *.intersection.bed 
rm "$i".path.*
#*****

done
