#!/bin/bash

for i in H3K4me1 H3K4me2 H3K4me3 H3K27ac H3K9ac H3K9me3 H3K27me3 H3K36me3 H4K20me1

do

cd ../analysis/all.marks/"$i"/

mkdir QN.merged

cd QN.merged

#*****
# 1. merge rep.1 and rep.2
paste <(tail -n+2 ../"$i".R1.matrix.tsv) <(tail -n+2 ../"$i".R2.matrix.tsv | cut -f2-) | sed '1igene_id\tH000.2T\tH003.2T\tH006.2T\tH009.2T\tH012.2T\tH018.2T\tH024.2T\tH036.2T\tH048.2T\tH072.2T\tH120.2T\tH168.2T\tH000.3T\tH003.3T\tH006.3T\tH009.3T\tH012.3T\tH018.3T\tH024.3T\tH036.3T\tH048.3T\tH072.3T\tH120.3T\tH168.3T' > "$i".R1.R2.matrix.tsv &&
#*****

#*****
# 2. run QN on merged matrix
Rscript ../../../bin/quantile.normalization.R "$i".R1.R2.matrix.tsv > "$i".R1.R2.matrix.after.QN.merged.tsv &&
#*****

#*****
# 3. split joint matrix into rep. 1 and rep. 2 matrices
awk 'BEGIN{FS=OFS="\t"}{if (NR==1){print "gene_id", $0} else {print $0}}' "$i".R1.R2.matrix.after.QN.merged.tsv | cut -f-13 | sed '1s/gene_id\t//' > "$i".R1.matrix.after.QN.merged.tsv &&
awk 'BEGIN{FS=OFS="\t"}{if (NR==1){print "gene_id", $0} else {print $0}}' "$i".R1.R2.matrix.after.QN.merged.tsv | cut -f1,14- | sed '1s/gene_id\t//' > "$i".R2.matrix.after.QN.merged.tsv &&
#*****

#*****
# 4. take the average of replicate 1 and 2
Rscript ../../../bin/matrix_matrix_mean.R -a "$i".R1.matrix.after.QN.merged.tsv -b "$i".R2.matrix.after.QN.merged.tsv -o "$i".matrix.after.QN.merged.tsv -r 2 &&
awk 'BEGIN{FS=OFS="\t"}{if (NR==1){for (i=1; i<=NF; i++){split($i, a, "."); $i=a[1]}; print $0} else {split($1, b, "."); $1=b[1]; print $0}}' "$i".matrix.after.QN.merged.tsv > tmp; mv tmp "$i".matrix.after.QN.merged.tsv &&
#*****

#*****
# 5. run maSigPro
Rscript ../../../bin/maSigPro.wrapper.R --rep1 "$i".R1.matrix.after.QN.merged.tsv --rep2 "$i".R2.matrix.after.QN.merged.tsv -o "$i".QN.merged.maSigPro.out.tsv &&
#*****

#*****
# 6. final histone mark matrix
/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/bin/selectMatrixRows.sh <(tail -n+2 ../../QN.merged/expression.QN.merged.maSigPro.out.tsv | cut -f1) "$i".matrix.after.QN.merged.tsv > "$i".matrix.tsv &&
#*****

#*****
# 7. number of genes that are variable out of the 8129 genes selected according to expression
grep -Ff <(tail -n+2 "$i".QN.merged.maSigPro.out.tsv | cut -f1 | cut -f1 -d ".") "$i".matrix.tsv | wc -l
#*****

done
