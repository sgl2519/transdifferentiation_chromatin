cd ../analysis/enhancers/ENCODE-rE2G.blood/

#!/bin/bash

#******
# 1. prepare bed file of coordinates
zcat interacting_pairs_filtered_100bp.merged_closest.tsv.gz | awk 'BEGIN{FS=OFS="\t"}{print $5, $6, $7, $4$5":"$6"-"$7, 0, ".", $4$5":"$6"-"$7}' > interacting_pairs_filtered_100bp.merged_closest.bed
#******

for i in H3K4me1 H3K4me2 H3K4me3 H3K27ac H3K9ac H3K9me3 H3K27me3 H4K20me1 H3K36me3

do

#******
# 2. prepare index file
mkdir all.marks/"$i"/
cd all.marks/"$i"/

ln -s ../../../analysis/all.marks/"$i"/"$i".R1.pvalue.bw.txt "$i".R1.pvalue.bw.txt
ln -s ../../../analysis/all.marks/"$i"/"$i".R2.pvalue.bw.txt "$i".R2.pvalue.bw.txt

ln -s ../../../analysis/all.marks/"$i"/"$i".peaks.txt "$i".peaks.txt 
#******

#******
# 6. prepare folder for bwtool analysis
mkdir bwtool 
#******

#******
# 3. define bedfile for this mark
myBed="../../../analysis/enhancers/ENCODE-rE2G.blood/interacting_pairs_filtered_100bp.merged_closest.bed" 
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

#conda deactivate

done

