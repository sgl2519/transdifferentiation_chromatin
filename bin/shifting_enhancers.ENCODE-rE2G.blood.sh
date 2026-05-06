#!/bin/bash

#*****
# Shift distal enhancer coordinates
#*****

cd ../analysis/enhancers/ENCODE-rE2G.blood
mkdir shifting
cd shifting

# 5 Kb shift
mkdir 5Kb
../../../../bin/join.py -b <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "chr"); print a[1], $0}' ../interacting_pairs_filtered_100bp.merged_closest.bed) -a <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "."); print a[1], $0}' ../../../..//bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed) | cut -f2- | awk 'BEGIN{FS=OFS="\t"}{if ($13=="+"){print $1, $2, $3, $4, $5, $6, $7, $9} else {print $1, $2, $3, $4, $5, $6, $7, $10}}' | awk 'BEGIN{FS=OFS="\t"}{if ($3<$8){$2=($2 - 5000); $3=($3 - 5000)} else {$2=($2 + 5000); $3=($3 + 5000)}; print $1, $2, $3, $4, $5, $6, $7}' > 5Kb/interacting_pairs_filtered_100bp.merged_closest.bed

# 10 Kb shift
mkdir 10Kb
../../../../bin/join.py -b <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "chr"); print a[1], $0}' ../interacting_pairs_filtered_100bp.merged_closest.bed) -a <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "."); print a[1], $0}' ../../../..//bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed) | cut -f2- | awk 'BEGIN{FS=OFS="\t"}{if ($13=="+"){print $1, $2, $3, $4, $5, $6, $7, $9} else {print $1, $2, $3, $4, $5, $6, $7, $10}}' | awk 'BEGIN{FS=OFS="\t"}{if ($3<$8){$2=($2 - 10000); $3=($3 - 10000)} else {$2=($2 + 10000); $3=($3 + 10000)}; print $1, $2, $3, $4, $5, $6, $7}' > 10Kb/interacting_pairs_filtered_100bp.merged_closest.bed

# 20 Kb shift
mkdir 20Kb
../../../../bin/join.py -b <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "chr"); print a[1], $0}' ../interacting_pairs_filtered_100bp.merged_closest.bed) -a <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "."); print a[1], $0}' ../../../..//bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed) | cut -f2- | awk 'BEGIN{FS=OFS="\t"}{if ($13=="+"){print $1, $2, $3, $4, $5, $6, $7, $9} else {print $1, $2, $3, $4, $5, $6, $7, $10}}' | awk 'BEGIN{FS=OFS="\t"}{if ($3<$8){$2=($2 - 20000); $3=($3 - 20000)} else {$2=($2 + 20000); $3=($3 + 20000)}; print $1, $2, $3, $4, $5, $6, $7}' > 20Kb/interacting_pairs_filtered_100bp.merged_closest.bed

# 50 Kb shift
mkdir 50Kb
../../../../bin/join.py -b <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "chr"); print a[1], $0}' ../interacting_pairs_filtered_100bp.merged_closest.bed) -a <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "."); print a[1], $0}' ../../../..//bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed) | cut -f2- | awk 'BEGIN{FS=OFS="\t"}{if ($13=="+"){print $1, $2, $3, $4, $5, $6, $7, $9} else {print $1, $2, $3, $4, $5, $6, $7, $10}}' | awk 'BEGIN{FS=OFS="\t"}{if ($3<$8){$2=($2 - 50000); $3=($3 - 50000)} else {$2=($2 + 50000); $3=($3 + 50000)}; print $1, $2, $3, $4, $5, $6, $7}' > 50Kb/interacting_pairs_filtered_100bp.merged_closest.bed
awk '$2>=0' 50Kb/interacting_pairs_filtered_100bp.merged_closest.bed > tmp; mv tmp 50Kb/interacting_pairs_filtered_100bp.merged_closest.bed

# 100 Kb shift
mkdir 100Kb
../../../../bin/join.py -b <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "chr"); print a[1], $0}' ../interacting_pairs_filtered_100bp.merged_closest.bed) -a <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "."); print a[1], $0}' ../../../..//bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed) | cut -f2- | awk 'BEGIN{FS=OFS="\t"}{if ($13=="+"){print $1, $2, $3, $4, $5, $6, $7, $9} else {print $1, $2, $3, $4, $5, $6, $7, $10}}' | awk 'BEGIN{FS=OFS="\t"}{if ($3<$8){$2=($2 - 100000); $3=($3 - 100000)} else {$2=($2 + 100000); $3=($3 + 100000)}; print $1, $2, $3, $4, $5, $6, $7}' > 100Kb/interacting_pairs_filtered_100bp.merged_closest.bed
awk '$2>=0' 100Kb/interacting_pairs_filtered_100bp.merged_closest.bed > tmp; mv tmp 100Kb/interacting_pairs_filtered_100bp.merged_closest.bed

# 1 Mb shift
mkdir 1Mb
../../../../bin/join.py -b <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "chr"); print a[1], $0}' ../interacting_pairs_filtered_100bp.merged_closest.bed) -a <(awk 'BEGIN{FS=OFS="\t"}{split($4, a, "."); print a[1], $0}' ../../../..//bed.files/gencode.v24.protein.coding.non.overlapping.genes.gene.body.bed) | cut -f2- | awk 'BEGIN{FS=OFS="\t"}{if ($13=="+"){print $1, $2, $3, $4, $5, $6, $7, $9} else {print $1, $2, $3, $4, $5, $6, $7, $10}}' | awk 'BEGIN{FS=OFS="\t"}{if ($3<$8){$2=($2 - 1000000); $3=($3 - 1000000)} else {$2=($2 + 1000000); $3=($3 + 1000000)}; print $1, $2, $3, $4, $5, $6, $7}' > 1Mb/interacting_pairs_filtered_100bp.merged_closest.bed
awk '$2>=0' 1Mb/interacting_pairs_filtered_100bp.merged_closest.bed > tmp; mv tmp 1Mb/interacting_pairs_filtered_100bp.merged_closest.bed
