#************
# LIBRARIES *
#************

library(ggplot2)
library(reshape2)

#********
# BEGIN *
#********
#*
# Plot activated and already-expressed genes ----

# 1. set working directory
setwd("../analysis/all.marks")

# 2. the marks we're analyzing
marks <- c("H3K27ac", "H3K9ac", "H4K20me1", "H3K4me3", "H3K4me1",
           "H3K36me3", "H3K4me2", "H3K9me3", "H3K27me3")

# 3. read marks' 5 decision-tree classes dataframes

## 3.1. read info of first mark (H3K27ac)
tmp <- read.table(paste0(marks[1], "/QN.merged/", marks[1], ".6.groups.tsv"), h=T, sep="\t", stringsAsFactors = F)

## 3.2. merge peak_not_TSS and no_peak into unmarked class
tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group)

x <- tmp[, "pearson_estimate_mark_level", drop = F]
colnames(x)[ncol(x)] <- marks[1]

## 3.3. add info regarding all other marks
for ( i in 2:9 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".6.groups.tsv"), h=T, sep="\t", stringsAsFactors = F)
  tmp$group <- gsub("peak_not_TSS", "no_peak", tmp$group)
  
  stopifnot(identical(rownames(x), rownames(tmp)))
  x <- cbind(x, tmp[, "pearson_estimate_mark_level"])
  colnames(x)[ncol(x)] <- marks[i]
  
}

rm(tmp)

# 4. genes not expressed at 0 hours
genes.257 <- read.table("expression/257.notExpressed.0h.txt",
                        stringsAsFactors = F)
colnames(genes.257) <- "gene_id"

# 5. subset cc df for 257 genes
x.257 <- x[rownames(x) %in% genes.257$gene_id, ]
x.257$type <- "257 activated genes"

# 6. genes expressed at 0 hours
genes.629 <- read.table("expression/629.expressed.0h.txt",
                        stringsAsFactors = F)
colnames(genes.629) <- "gene_id"

# 7. subset cc df for 629 genes
x.629 <- x[rownames(x) %in% genes.629$gene_id, ]
x.629$type <- "629 already expr. genes"

x <- rbind(x.257, x.629)
x.melt <- reshape2::melt(x)

x.melt$variable <- factor(x.melt$variable, levels = c("H3K4me1",
                                                      "H3K4me2",
                                                      "H3K27ac",
                                                      "H3K9ac",
                                                      "H3K4me3",
                                                      "H3K36me3",
                                                      "H4K20me1",
                                                      "H3K27me3",
                                                      "H3K9me3"
))

palette <- c("H3K9ac" = "#AF4D85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#D199B9",
             "H3K27me3" = "#1D2976",
             "H3K9me3" = "#A7ADD4",
             "H3K36me3" = "#7FBC41",
             "H4K20me1" = "#4C7027",
             "H3K4me1" = "#E5AB00",
             "H3K4me2" = "#A67C00")

pdf("fig_4a.pdf",
    width = 6, height = 3)
ggplot(x.melt, aes(x=variable, y=value, fill=variable)) +
  geom_violin(alpha=.4, color = "white", scale='width') +
  geom_boxplot(width=.25, alpha = .8, outlier.shape = NA) +
  # geom_boxplot(alpha = .7) + #, outlier.shape = NA, width=0.5) +
  # geom_violin(alpha=.7, colour="white", width = 1) +
  facet_grid(.~type) +
  scale_fill_manual(values = palette) +
  guides(fill='none') +
  ylab(expression(paste("Pearson ", italic("r")))) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        # axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=13),
        axis.text.x = element_text(size=13, hjust=1, vjust=.5, angle=90),
        panel.border = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # plot.title = element_text(size = 15, hjust = .5),
        plot.title = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 15, angle = -360),
        panel.spacing = unit(0.1, "lines")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
  labs(title = "promoters/gene bodies")

dev.off()
