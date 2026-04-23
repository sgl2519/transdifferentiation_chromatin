#!/bin/bash
#*******
# 1. generate folders to store results
mkdir out
#*******

#*******
# 2. retrieve list of
# files with ENCODE peak calls
ls path/to/peak_calls/*bed.gz | sort -k3,3 -k1,1 > peaks.files.txt
#*******

#*******
# 3. generate folder to store files
# of selected genes for each mark
mkdir selected.genes
#*******

#*******
# 4. retrieve selected genes for each mark
# aka genes with peaks of the mark 
# in (+- 5Kb promoter region) in all 12 time-points
for mark in H3K4me3 H3K4me2 H3K4me1 H3K36me3 H4K20me1 H3K9ac H3K9me3 H3K27me3 H3K27ac; do \
grep "$mark" peaks.files.txt | while read file; do bedtools intersect -a ../bed.files/gencode.v24.protein.coding.non.overlapping.genes.5Kb.upstream.downstream.TSS.bed -b $file -u; done \
| sort | uniq -c | awk '$1==12{print $5}' > selected.genes/"$mark".selected.genes.txt; \
done
#*******

#*******
# 5. prepare folder to store
# paths to bw files for each mark and replicate
mkdir bw.files
#*******

#*******
# 6. retrieve paths to bw files 
# for each mark and replicate
for mark in H3K4me3 H3K4me2 H3K4me1 H3K36me3 H4K20me1 H3K9ac H3K9me3 H3K27me3 H3K27ac
do 
for rep in 1 2
do 
grep -F "$mark""X""$rep" ../references/pipeline.db | sort -k1,1 | cut -f2 > bw.files/"$mark".R"$rep".pvalue.bw.txt
done
done
#******

#******
# 7. aggregation plots - tss (5 Kb up/down-stream)
for mark in H3K9ac H3K27ac H3K4me1 H3K4me2 H3K36me3 H3K9me3 H4K20me1 H3K27me3 H3K4me3
do
# 7.1. prepare tmp bed file
grep -Ff selected.genes/"$mark".selected.genes.txt ../bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed > tmp.bed

# 7.2. get aggregate tables
for rep in R1 R2
do
../bin/bwtool.aggregate.ChIPseq.sh --bw bw.files/"$mark"."$rep".pvalue.bw.txt --bedfile tmp.bed --type tss --outFile out/"$mark".aggregation.plot."$rep".tsv
done
        
# 7.3. rm tmp bed file
rm tmp.bed
done
#******

#******
# 8. aggregation plots - gene body (5 Kb up/down-stream)
for mark in H3K36me3 H4K20me1
do
# 8.1. prepare tmp bed file
grep -Ff selected.genes/"$mark".selected.genes.txt ../bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed > tmp.bed

# 8.2. retrieve meta value
pos=$(awk 'function ceil(x, y){y=int(x); return(x>y?y+1:y)}{n++}END{$0=n/2; print ceil($0)}' tmp.bed)
myMeta=$(awk '{print $3-$2}' tmp.bed | sort -n | awk -v pos="$pos" 'NR==pos{print}')

# 8.3. get aggregate tables
for rep in R1 R2
do
../bin/bwtool.aggregate.ChIPseq.sh --bw bw.files/"$mark"."$rep".pvalue.bw.txt --bedfile tmp.bed --type gene_body --meta "$myMeta" --outFile out/"$mark".aggregation.plot."$rep".gene.body.tsv
done
        
# 8.4 rm tmp bed file
rm tmp.bed
done
#******
