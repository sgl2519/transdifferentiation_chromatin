#************
# LIBRARIES *
#************

library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(circlize)

# Plot the normalized signal of marking ----

#************
# FUNCTIONS *
#************

function1 <- function(mark) {
  
  # 1. read mark matrix
  mark.matrix <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.tsv"), h=T, sep="\t")
  
  # 2. reorder mark matrix according to metadata
  mark.matrix <- mark.matrix[rownames(metadata), ]
  stopifnot(identical(rownames(metadata), rownames(mark.matrix)))
  
  # 3. center and scale mark matrix
  mark.matrix.scaled <- as.data.frame(t(scale(t(mark.matrix))))
  colnames(mark.matrix.scaled) <- colnames(expression.matrix.scaled)
  
  # 4. read matrix of 6 groups
  # and merge peak_not_TSS and no_peak
  mark.6.groups <- read.table(paste0(mark, "/QN.merged/", mark, ".6.groups.tsv"),
                              h=T, sep="\t")
  mark.6.groups$group <- gsub("peak_not_TSS", "no_peak", mark.6.groups$group)
  
  
  # 5. genes w/o peaks of the mark or stable
  mark.no.peak.ns.genes <- rownames(mark.6.groups[mark.6.groups$group %in% c("no_peak", "stable"), ])
  mark.metadata.no.peak.ns.genes <- metadata[rownames(metadata) %in% mark.no.peak.ns.genes, ]
  stopifnot(identical(rownames(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.no.peak.ns.genes), ]),
                      rownames(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.no.peak.ns.genes), ])))
  
  
  # 6. genes w/ peak of the mark and variable profile
  mark.sig.genes <- rownames(mark.6.groups[mark.6.groups$group %in% c("positively_correlated",
                                                                      "negatively_correlated",
                                                                      "not_correlated"), ])
  mark.metadata.sig.genes <- metadata[rownames(metadata) %in% mark.sig.genes, ]
  stopifnot(identical(rownames(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ]),
                      rownames(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ])))
  
  mark.partition1 = rbind(mark.metadata.sig.genes,
                          mark.metadata.no.peak.ns.genes)[, "final_class"]
  
  
  mark.ht_list = Heatmap(mark.partition1, 
                         name = "profiles", 
                         col= palette[names(palette) %in% c("bending",
                                                            "down-regulated",
                                                            "peaking",
                                                            "up-regulated")],
                         show_heatmap_legend = F,
                         show_row_names = FALSE, width = unit(8, "mm"),
                         show_column_names = F,
                         split = c(rep("D", nrow(mark.metadata.sig.genes)), 
                                   rep("S", nrow(mark.metadata.no.peak.ns.genes))),
                         gap = unit(3, "mm"),
                         row_title_gp = gpar(fontsize = 26)) +
    
    Heatmap(rbind(expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ], 
                  expression.matrix.scaled[rownames(expression.matrix.scaled) %in% rownames(mark.metadata.no.peak.ns.genes), ]),
            col = palette2,
            name = "z-score", column_title = "expression",
            column_title_gp = gpar(fontsize = 29),
            column_names_gp = gpar(fontsize = 46),
            show_row_names = FALSE, width = unit(50, "mm"),
            cluster_rows = F, cluster_columns = F,
            show_column_names = F,
            show_heatmap_legend = F,
            heatmap_legend_param = list(title = "z-score"),
            split = c(rep("D", nrow(mark.metadata.sig.genes)), 
                      rep("S", nrow(mark.metadata.no.peak.ns.genes))),
            gap = unit(3, "mm"),
            border = T) +
    
    Heatmap(rbind(mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.sig.genes), ], 
                  mark.matrix.scaled[rownames(mark.matrix.scaled) %in% rownames(mark.metadata.no.peak.ns.genes), ]),
            col = palette2, 
            name = mark, column_title = mark,
            column_title_gp = gpar(fontsize = 29),
            column_names_gp = gpar(fontsize = 46),
            show_row_names = FALSE, width = unit(50, "mm"),
            cluster_rows = F, cluster_columns = F,
            show_column_names = F,
            show_heatmap_legend = F,
            heatmap_legend_param = list(title = "z-score"),
            split = c(rep("D", nrow(mark.metadata.sig.genes)), 
                      rep("U", nrow(mark.metadata.no.peak.ns.genes))),
            gap = unit(3, "mm"),
            border = T)
  
  return(mark.ht_list)

  
}


#********
# BEGIN *
#********

setwd("../analysis/all.marks/")

palette <- c("down-regulated" = "#2d7f89", 
             "bending" = "#7acbd5",
             "up-regulated" = "#89372d",
             "peaking" = "#d5847a"
)

palette2 = rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))


# 1. import expression matrix
expression.matrix <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")

# 2. import metadata and check row order is the same as in expression matrix
metadata <- read.table("expression/QN.merged/metadata.class2.tsv", h=T, sep="\t")
stopifnot(identical(rownames(expression.matrix), rownames(metadata)))

# 3. work with final_class column
metadata$final_class <- gsub("regulation", "-regulated", metadata$final_class)

# 4. add to metadata info about average expression and time-point
metadata$avg_exp <- apply(expression.matrix, 1, mean)
metadata$tp1 <- apply(expression.matrix, 1, which.max)
metadata$tp2 <- apply(expression.matrix, 1, which.min)
metadata$tp <- ifelse(metadata$final_class %in% c("bending", "down-regulated"), 
                      metadata$tp2, 
                      metadata$tp1)
metadata <- metadata[order(metadata$final_class, metadata$tp, metadata$avg_exp), ]

# 5. reorder expression matrix according to metadata row order
expression.matrix <- expression.matrix[rownames(metadata), ]
stopifnot(identical(rownames(expression.matrix), rownames(metadata)))

# 6. center and scale expression matrix
expression.matrix.scaled <- as.data.frame(t(scale(t(expression.matrix))))
colnames(expression.matrix.scaled) <- c("0", "3", "6", "9", "12",
                                        "18", "24", "36", "48", "72", "120", "168")

# 7. the marks we're analyzing
marks <- c("H3K9ac", "H3K27ac", "H4K20me1", "H3K4me3", "H3K4me1", "H3K4me2", "H3K36me3", "H3K27me3", "H3K9me3")


# 8. make plots
lop <- list()

for ( i in 1:9 ) {
  
  lop[[i]] <- function1(mark = marks[i])
  
}

# 9. extract legend
foo1 = color_mapping_legend(lop[[9]]@ht_list[[1]]@matrix_color_mapping, plot = FALSE) 
foo2 = color_mapping_legend(lop[[9]]@ht_list[[2]]@matrix_color_mapping, plot = FALSE) 


# 9. arrange heatmaps in grid

pdf("../../figures/fig_S7b.pdf",
    width = 28, height = 17)

grid.newpage()
pushViewport(viewport(layout = grid.layout(nr = 2, nc = 5)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(lop[[1]], newpage = FALSE, padding = unit(c(10, 10, 10, 10), "mm"))
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(lop[[2]], newpage = FALSE, padding = unit(c(10, 10, 10, 10), "mm"))
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
draw(lop[[3]], newpage = FALSE, padding = unit(c(10, 10, 10, 10), "mm"))
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 4))
draw(lop[[4]], newpage = FALSE, padding = unit(c(10, 10, 10, 10), "mm"))
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 5))
draw(lop[[5]], newpage = FALSE, padding = unit(c(10, 10, 10, 10), "mm"))
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
draw(lop[[6]], newpage = FALSE, padding = unit(c(10, 10, 10, 10), "mm"))
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 2))
draw(lop[[7]], newpage = FALSE, padding = unit(c(10, 10, 10, 10), "mm"))
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 3))
draw(lop[[8]], newpage = FALSE, padding = unit(c(10, 10, 10, 10), "mm"))
upViewport()

pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 4))
draw(lop[[9]], newpage = FALSE, padding = unit(c(10, 10, 10, 10), "mm"))
upViewport()


# add legend
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 5))
draw(foo1, just = "top")
draw(foo2, just = "bottom")
upViewport()


upViewport()

dev.off()
