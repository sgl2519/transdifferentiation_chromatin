#************
# LIBRARIES *
#************

library(seriation)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(tidyr)

# Obtain sets of H3K27me3 pre-marked and not pre-marked genes
x <- read.table("expression/257.notExpressed.0h.txt",
                stringsAsFactors = F)

# Load matrix of H3K27me3 marking
h3k27me3 <- read.table('H3K27me3/QN.merged/ant.del.analysis/H3K27me3.peaks.dynamics.tsv')
h3k27me3$V1 <- gsub('\\..*','', h3k27me3$V1)
h3k27me3$V2 <- gsub(';.*','', h3k27me3$V2)

act.genes.premarked <- x[x$V1 %in% h3k27me3[h3k27me3$V2 == "H000" |
                                              h3k27me3$V2 == "H003", ]$V1, ]
act.genes.unmarked <- x[x$V1 %in% setdiff(x$V1, act.genes.premarked), ]

# 1. set working directory
setwd("../analysis/all.marks")

# 2. read set of 257 activated genes
x <- read.table("expression/257.notExpressed.0h.txt",
                stringsAsFactors = F)
colnames(x) <- "gene_id"

# 3. the histone marks we analyze
marks <- c("H3K4me1", "H3K4me2", "H3K27ac",
           "H3K9ac","H3K4me3", "H3K36me3",
           "H4K20me1", "H3K9me3", "H3K27me3")

# 4. list of time-points
hours <- c("H000", "H003", "H006", "H009","H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")

# 5. for every mark, read dataframe with info
# of presence/absence of peaks across time-points
# for the 257 activated genes
for ( i in 1:9) {
  
  # 5.1. read dataframe with
  # 1st coln = gene_id
  # 2nd coln = time-points with peaks
  
  # Read peak presence/absence matrices
  tmp <- read.table(paste0(marks[i],
                           "/QN.merged/ant.del.analysis/",
                           marks[i],
                           ".peaks.dynamics.tsv"),
                    stringsAsFactors = F)
  
  colnames(tmp) <- c("gene_id", "tp")
  
  # 5.2. fix gene_id (aka remove everything after dot)
  tmp$gene_id <- gsub("\\..*","", tmp$gene_id)
  
  # 5.3. create empty dataframe
  m <- data.frame(stringsAsFactors = F)
  
  # 5.4. for each of the 257 activated genes
  for (gene in x$gene_id) {
    
    # 5.4.1. create a vector of 0s (n=12)
    v <- rep(0, 12)
    if (gene %in% tmp$gene_id) {
      
      # 5.4.2. retrieve the time-points in which the gene has a peak
      stn <- as.character(unlist(strsplit(tmp[tmp$gene_id == gene, "tp"], ";")))
      
      # 5.4.3. assign v=1 to the time-points in which the gene has a peak
      v[which(hours %in% stn)] <- 1
      
    }
    m <- rbind(m, v)
    
  }
  
  peaks <- m
  
  # Read signal data
  tmp <- read.table(paste0(marks[i],
                           "/QN.merged/",
                           marks[i],
                           ".matrix.after.QN.merged.tsv"),
                    stringsAsFactors = F)
  
  # 5.3. create empty dataframe
  m <- tmp[x$gene_id, ]
  
  colnames(m) <- paste(hours, marks[i], sep="_")
  
  # Normalize only the signal on cases where there is peak
  for ( j in 1:nrow(peaks) ) {
    print(paste0(j, ': ', length(peaks[j, peaks[j, ] == 1])))
    
    if ( length(peaks[j, peaks[j, ] == 1]) > 1 ) {
      m[j, peaks[j, ] == 1] <- (m[j, peaks[j, ] == 1] - min(m[j, peaks[j, ] == 1]))/max(m[j, peaks[j, ] == 1])
      m[j, peaks[j, ] == 0] <- 0
      
    } else if ( length(peaks[j, peaks[j, ] == 1]) == 1 ) {
      m[j, ] <- rep(0, 12)
      
    } else {
      m[j, ] <- rep(0, 12)
      
    }
    
  }
  
  # # Normalize dividing subtracting the minimum and dividing by the maximum
  # m <- t(apply(m, 1, function(a) (a-min(a))/max(a) ))
  
  x <- cbind(x, m)
  
}


# 6. add rownames and remove gene_id column
rownames(x) <- x$gene_id
x$gene_id <- NULL


# 7. repeat the same procedure for expression
# setting 1 = tp at which expression > 1 TPM
expression.matrix <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
expression.matrix <- expression.matrix[rownames(expression.matrix) %in% rownames(x), ]

# 5.3. create empty dataframe
m <- expression.matrix[rownames(x), ]

# # Normalize
# m <- as.data.frame(t(scale(t(m))))

# Normalize dividing subtracting the minimum and dividing by the maximum
m <- t(apply(m, 1, function(a) (a-min(a))/max(a) ))

colnames(m) <- paste(hours, "expression", sep="_")
x <- cbind(x, m)


# 8. reorder columns in x according to order of marks
x <- x[, paste(rep(hours, 10), rep(c("H3K4me1", "H3K4me2",
                                     "expression",
                                     "H3K27ac", "H3K9ac",
                                     "H3K4me3",
                                     "H3K36me3", "H4K20me1",
                                     "H3K9me3", "H3K27me3"), each = 12), sep = "_")]

# Order based on time point of gene activation separated by H3K27me3 pre-marked and not pre-marked independently
## Extract the time point of gene activation for each gene
tp_activation <- apply(expression.matrix, 1, function(a) which(a > 1)[1])

ordered <- data.frame()
for ( names in list(act.genes.unmarked, act.genes.premarked) ) {
  
  tmp <- x[names, ]
  
  # Loop through the different time points of gene activation
  for ( t in 2:11 ) {
    if ( length(intersect(names, names(tp_activation[tp_activation == t]))) > 1 ) {
      tmp_ss <- tmp[intersect(names, names(tp_activation[tp_activation == t])), ]
      ser.obj <- seriate(as.matrix(tmp_ss[, 25:36]), method = "PCA_angle")
      x.ser <- tmp_ss[get_order(ser.obj, 1), ]
      
      ordered <- rbind(ordered, x.ser)
    } else {
      ordered <- rbind(ordered, tmp[intersect(names, names(tp_activation[tp_activation == t])), ])
      
    }
    
  }
  
}

x.ser <- x.ser[nrow(x.ser):1, ]
pdf("fig_4e.pdf", 
   width = 25, height = 7)
Heatmap(ordered[, 1:12],
        # col = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
        #             '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),
        # col = rev(c('#E5AB00', 'white')),
        # col = c("#fee6ce", "#ff1919"),
        col = colorRamp2(c(0, 0.001, 1), c("white", "#fee6ce", "#ff1919")),
        name = "H3K4me1", column_title = "H3K4me1",
        column_title_gp = gpar(fontsize = 29),
        column_names_gp = gpar(fontsize = 46),
        show_row_names = FALSE, width = unit(50, "mm"),
        cluster_rows = F, cluster_columns = F,
        show_column_names = F,
        show_heatmap_legend = F,
        heatmap_legend_param = list(title = "z-score"),
        row_split = c(rep('a_unmarked', length(act.genes.unmarked)),
                      rep('b_premarked', length(act.genes.premarked))),
        # split = c(rep("a", 165), 
        #           rep("b", 92)),
        gap = unit(3, "mm"),
        use_raster = TRUE,
        border = T) +
  
  Heatmap(ordered[, 13:24],
          # col = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
          #             '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),
          # col = rev(c('#E5AB00', 'white')),
          # col = c("#fee6ce", "#ff1919"),
          col = colorRamp2(c(0, 0.001, 1), c("white", "#fee6ce", "#ff1919")),
          name = "H3K4me2", column_title = "H3K4me2",
          column_title_gp = gpar(fontsize = 29),
          column_names_gp = gpar(fontsize = 46),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          row_split = c(rep('a_unmarked', length(act.genes.unmarked)),
                        rep('b_premarked', length(act.genes.premarked))),
          # split = c(rep("a", 165), 
          #           rep("b", 92)),
          gap = unit(3, "mm"),
          use_raster = TRUE,
          border = T) +
  
  Heatmap(ordered[, 37:48],
          # col = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
          #             '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),
          # col = rev(c('#E5AB00', 'white')),
          # col = c("#fee6ce", "#ff1919"),
          col = colorRamp2(c(0, 0.001, 1), c("white", "#fee6ce", "#ff1919")),
          name = "H3K27ac", column_title = "H3K27ac",
          column_title_gp = gpar(fontsize = 29),
          column_names_gp = gpar(fontsize = 46),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          row_split = c(rep('a_unmarked', length(act.genes.unmarked)),
                        rep('b_premarked', length(act.genes.premarked))),
          # split = c(rep("a", 165), 
          #           rep("b", 92)),
          gap = unit(3, "mm"),
          use_raster = TRUE,
          border = T) +
  
  Heatmap(ordered[, 49:60],
          # col = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
          #             '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),
          # col = rev(c('#E5AB00', 'white')),
          # col = c("#fee6ce", "#ff1919"),
          col = colorRamp2(c(0, 0.001, 1), c("white", "#fee6ce", "#ff1919")),
          name = "H3K9ac", column_title = "H3K9ac",
          column_title_gp = gpar(fontsize = 29),
          column_names_gp = gpar(fontsize = 46),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          row_split = c(rep('a_unmarked', length(act.genes.unmarked)),
                        rep('b_premarked', length(act.genes.premarked))),
          # split = c(rep("a", 165), 
          #           rep("b", 92)),
          gap = unit(3, "mm"),
          use_raster = TRUE,
          border = T) +
  
  Heatmap(ordered[, 25:36],
          # col = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
          #             '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),
          # col = rev(c('#E5AB00', 'white')),
          # col = c("#ddfafd", "#3182bd"),
          col = colorRamp2(c(0, 0.001, 1), c("white", "#ddfafd", "#3182bd")),
          name = "expression", column_title = "expression",
          column_title_gp = gpar(fontsize = 29),
          column_names_gp = gpar(fontsize = 46),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          row_split = c(rep('a_unmarked', length(act.genes.unmarked)),
                        rep('b_premarked', length(act.genes.premarked))),
          # split = c(rep("a", 165),
          #           rep("b", 92)),
          gap = unit(3, "mm"),
          use_raster = TRUE,
          border = T) +
  
  Heatmap(ordered[, 61:72],
          # col = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
          #             '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),
          # col = rev(c('#E5AB00', 'white')),
          # col = c("#fee6ce", "#ff1919"),
          col = colorRamp2(c(0, 0.001, 1), c("white", "#fee6ce", "#ff1919")),
          name = "H3K4me3", column_title = "H3K4me3",
          column_title_gp = gpar(fontsize = 29),
          column_names_gp = gpar(fontsize = 46),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          row_split = c(rep('a_unmarked', length(act.genes.unmarked)),
                        rep('b_premarked', length(act.genes.premarked))),
          # split = c(rep("a", 165), 
          #           rep("b", 92)),
          gap = unit(3, "mm"),
          use_raster = TRUE,
          border = T) +
  
  Heatmap(ordered[, 73:84],
          # col = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
          #             '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),
          # col = rev(c('#E5AB00', 'white')),
          # col = c("#fee6ce", "#ff1919"),
          col = colorRamp2(c(0, 0.001, 1), c("white", "#fee6ce", "#ff1919")),
          name = "H3K36me3", column_title = "H3K36me3",
          column_title_gp = gpar(fontsize = 29),
          column_names_gp = gpar(fontsize = 46),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          row_split = c(rep('a_unmarked', length(act.genes.unmarked)),
                        rep('b_premarked', length(act.genes.premarked))),
          # split = c(rep("a", 165), 
          #           rep("b", 92)),
          gap = unit(3, "mm"),
          use_raster = TRUE,
          border = T) +
  
  Heatmap(ordered[, 85:96],
          # col = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
          #             '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),
          # col = rev(c('#E5AB00', 'white')),
          # col = c("#fee6ce", "#ff1919"),
          col = colorRamp2(c(0, 0.001, 1), c("white", "#fee6ce", "#ff1919")),
          name = "H4K20me1", column_title = "H4K20me1",
          column_title_gp = gpar(fontsize = 29),
          column_names_gp = gpar(fontsize = 46),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          row_split = c(rep('a_unmarked', length(act.genes.unmarked)),
                        rep('b_premarked', length(act.genes.premarked))),
          # split = c(rep("a", 165), 
          #           rep("b", 92)),
          gap = unit(3, "mm"),
          use_raster = TRUE,
          border = T) +
  
  Heatmap(ordered[, 97:108],
          # col = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
          #             '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),
          # col = rev(c('#E5AB00', 'white')),
          # col = c("#fee6ce", "#ff1919"),
          col = colorRamp2(c(0, 0.001, 1), c("white", "#fee6ce", "#ff1919")),
          name = "H3K9me3", column_title = "H3K9me3",
          column_title_gp = gpar(fontsize = 29),
          column_names_gp = gpar(fontsize = 46),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          row_split = c(rep('a_unmarked', length(act.genes.unmarked)),
                        rep('b_premarked', length(act.genes.premarked))),
          # split = c(rep("a", 165), 
          #           rep("b", 92)),
          gap = unit(3, "mm"),
          use_raster = TRUE,
          border = T) +
  
  Heatmap(ordered[, 109:120],
          # col = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
          #             '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')),
          # col = rev(c('#E5AB00', 'white')),
          # col = c("#fee6ce", "#ff1919"),
          col = colorRamp2(c(0, 0.001, 1), c("white", "#fee6ce", "#ff1919")),
          name = "H3K27me3", column_title = "H3K27me3",
          column_title_gp = gpar(fontsize = 29),
          column_names_gp = gpar(fontsize = 46),
          show_row_names = FALSE, width = unit(50, "mm"),
          cluster_rows = F, cluster_columns = F,
          show_column_names = F,
          show_heatmap_legend = F,
          heatmap_legend_param = list(title = "z-score"),
          row_split = c(rep('a_unmarked', length(act.genes.unmarked)),
                        rep('b_premarked', length(act.genes.premarked))),
          # split = c(rep("a", 165), 
          #           rep("b", 92)),
          gap = unit(3, "mm"),
          use_raster = TRUE,
          border = T)
dev.off()
