#************
# LIBRARIES *
#************

library(dplyr)
library(tidyr)
library(seriation)

#************
# FUNCTIONS *
#************

function1 <- function(mark) {

  # 1. read mark matrix of 12448 genes after QN merged normalization
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                            h=T, sep="\t")

  # 2. keep only stably expressed genes (aka genes in expression matrix)
  mark.matrix <- mark.matrix[rownames(mark.matrix) %in% rownames(expression.matrix.scaled), ]

  # 3. check order of rows is the same as with expression
  mark.matrix <- mark.matrix[rownames(expression.matrix), ]
  stopifnot(identical(rownames(expression.matrix.scaled), rownames(mark.matrix)))

  # 4. center and scale mark matrix of stably expressed genes
  mark.matrix.scaled <- as.data.frame(t(scale(t(mark.matrix))))

  # 5. retrieve set of marked genes
  marked.genes <- read.table(paste0(mark, "/QN.merged/all.genes.intersecting.peaks.tsv"),
                             h=F, sep="\t", stringsAsFactors = F)
  marked.genes <- marked.genes$V1

  # 6. retrieve set of genes w/ variable profile for the mark
  mark.sig.genes <- read.table(paste0(mark, "/QN.merged/", mark, ".QN.merged.maSigPro.out.tsv"),
                               h=T, sep="\t")
  mark.sig.genes$gene_id <- rownames(mark.sig.genes)
  mark.sig.genes <- mark.sig.genes %>% separate(gene_id, c("gene", "id"), "\\.")
  rownames(mark.sig.genes) <- mark.sig.genes$gene
  mark.sig.genes$gene <- NULL
  mark.sig.genes$id <- NULL
  mark.sig.genes <- rownames(mark.sig.genes)[rownames(mark.sig.genes) %in% rownames(expression.matrix.scaled)]

  # 7. retrieve differentially marked genes (aka w/ peak of the mark and variable profile)
  mark.sig.genes <- intersect(mark.sig.genes, marked.genes)

  print(mark)
  print(length(mark.sig.genes))
  print(length(setdiff(intersect(marked.genes, rownames(expression.matrix.scaled)), mark.sig.genes)))

  # 8. return mark matrix scaled for differentially marked genes
  x <- mark.matrix.scaled[rownames(mark.matrix.scaled) %in% mark.sig.genes, ]
  x$mark <- mark
  x$gene_id <- rownames(x)

  return(x)

}


#********
# BEGIN *
#********

# 1. set working directory
setwd("analysis/all.marks/")

# 2. read expression matrix of 12448 genes
expression.matrix <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                                h=T, sep="\t")

# 3. retrieve set of not expressed genes
not.expressed.genes <- read.table("expression/silent.genes.txt", h=F, sep="\t",
                                  stringsAsFactors = F)
not.expressed.genes <- not.expressed.genes$V1


# 4. retrieve set of DE genes
DE.genes <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
DE.genes <- rownames(DE.genes)


# 5. keep expression matrix only for stably expressed genes
expression.matrix <- expression.matrix[!(rownames(expression.matrix) %in%
                                           c(not.expressed.genes, DE.genes)), ]

# 6. center and scale expression matrix of flat genes
expression.matrix.scaled <- as.data.frame(t(scale(t(expression.matrix))))

# 7. the marks we're analyzing
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K36me3", "H3K4me3", "H3K4me1",
           "H3K4me2", "H3K9me3", "H3K27me3")

# 8. retrieve matrix of scaled and centered mark values for differentially marked genes
m <- data.frame(stringsAsFactors = F)

for ( i in 1:9 ) {
  m <- rbind(m, function1(mark = marks[i]))

}

# 9. apply seriation
ser.obj <- seriate(as.matrix(m[, 1:12]), method = "PCA_angle")
m.ser <- m[get_order(ser.obj, 1), ]

# 10. prepare df for heatmap
df.plot <- data.frame(stringsAsFactors = F)
y <- c()

for ( i in 1:9 ) {

  # 10.1. mark matrix reordered according to seriation order
  tmp1 <- m.ser[m.ser$mark == marks[i], ]
  rownames(tmp1) <- tmp1$gene_id # recover gene_id
  tmp1 <- tmp1[, 1:12]
  colnames(tmp1) <- paste(colnames(tmp1), "mark", sep="_")

  # 10.2. expression matrix reordered according to seriation order of mark
  tmp2 <- expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(tmp1), ]
  tmp2 <- tmp2[rownames(tmp1), ]
  stopifnot(identical(rownames(tmp1), rownames(tmp2)))
  colnames(tmp2) <- paste(colnames(tmp2), "expression", sep="_")

  # 10.3. store number of genes for each mark
  y <- c(y, nrow(tmp1))

  # 10.4. store df for mark and expression
  df.plot <- rbind(df.plot, cbind(tmp1, tmp2))

}

# 11. prepare mark partition for heatmap
mark.partition <- c()
for ( i in 1:9 ) {
  mark.partition <- c(mark.partition, rep(as.character(y[i]), y[i]))

}

mark.partition <- factor(mark.partition, levels = as.character(y))

# 12. define color palettes
palette1 <- c("#630039","#af4d85", "#4c7027",
              "#7fbc41", "#d199b9","#e5ab00",
              "#a67c00","#a7add4","#1d2976")

names(palette1) <- levels(mark.partition)

palette2 = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))

#************
# LIBRARIES *
#************

library(ComplexHeatmap)
library(circlize)

# 2. make heatmap

pdf("figures/fig_S7c.pdf", 
    height=18, width=13)

ht_list <- Heatmap(mark.partition, 
                   name = "marks", 
                   col= palette1,
                   show_heatmap_legend = F,
                   show_row_names = FALSE, width = unit(15, "mm"),
                   show_column_names = F,
                   split = mark.partition,
                   gap = unit(4, "mm"),
                   row_title_gp = gpar(fontsize = 45),
                   row_title_rot = 0) +
  
  Heatmap(df.plot[, 13:24],
          col = palette2,
          name = "log2(TPM+1)", column_title = "expression",
          column_title_gp = gpar(fontsize = 50),
          show_row_names = FALSE, width = unit(120, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = T,
          heatmap_legend_param = list(title = "z-score",
                                      title_gp = gpar(fontsize=45, fontface="bold"),
                                      labels_gp = gpar(fontsize=45),
                                      grid_width = unit(12, "mm"),
                                      grid_height = unit(12, "mm"),
                                      grid_border = unit(12, "mm"),
                                      legend_direction = "horizontal",
                                      title_position = "topcenter"),
          split = mark.partition,
          gap = unit(4, "mm"),
          border = T) +
  
  Heatmap(df.plot[, 1:12],
          col = palette2, 
          name = "mark", column_title = "histone mark",
          column_title_gp = gpar(fontsize = 50),
          show_row_names = FALSE, width = unit(120, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          split = mark.partition,
          gap = unit(4, "mm"),
          border = T)

draw(ht_list, padding = unit(c(10, 10, 10, 10), "mm"),
     heatmap_legend_side = "bottom")

dev.off()
