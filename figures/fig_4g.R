#************
# LIBRARIES *
#************

# library(seriation)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(tidyr)

# Promoters ----

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
      
      # 5.4.3. assign v=i only to the first time-point in which we have a peak
      # v[min(which(hours %in% stn)):12] <- 1
      v[min(which(hours %in% stn))] <- i
      
    }
    m <- rbind(m, v)
    
  }
  
  colnames(m) <- paste(hours, marks[i], sep="_")
  x <- cbind(x, m)
  
}


# 6. add rownames and remove gene_id column
rownames(x) <- x$gene_id
x$gene_id <- NULL


# # 7. repeat the same procedure for expression 
# setting 11 = tp at which expression > 1 TPM
expression.matrix <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
expression.matrix <- expression.matrix[rownames(expression.matrix) %in% rownames(x), ]

m <- data.frame(stringsAsFactors = F)

for (gene in rownames(x)) {
  
  v <- rep(0, 12)
  stn <- as.numeric(expression.matrix[gene, ])
  
  # assign v=11 only to the first time-point in which we have expression > 1
  v[min(which(stn > 1))] <- 11
  
  m <- rbind(m, v)
  
}


colnames(m) <- paste(hours, "expression", sep="_")
x <- cbind(x, m)




# 8. reorder columns in x according to order of marks
x <- x[, paste(rep(hours, 10), rep(c("H3K4me1", "H3K4me2",
                                     "expression",
                                     "H3K27ac", "H3K9ac", 
                                     "H3K4me3",
                                     "H3K36me3", "H4K20me1", 
                                     "H3K9me3", "H3K27me3"), each = 12), sep = "_")]


# 9. retrieve time-point of mark deposition / beginning of expression
x$H3K4me1 <- apply(x[, 1:12], 1, function(x){min(which(x>0))})
x$H3K4me2 <- apply(x[, 13:24], 1, function(x){min(which(x>0))})
x$expression_tp <- apply(x[, 25:36], 1, function(x){min(which(x>0))})
x$H3K27ac <- apply(x[, 37:48], 1, function(x){min(which(x>0))})
x$H3K9ac <- apply(x[, 49:60], 1, function(x){min(which(x>0))})
x$H3K4me3 <- apply(x[, 61:72], 1, function(x){min(which(x>0))})
x$H3K36me3 <- apply(x[, 73:84], 1, function(x){min(which(x>0))})
x$H4K20me1 <- apply(x[, 85:96], 1, function(x){min(which(x>0))})
x$H3K9me3 <- apply(x[, 97:108], 1, function(x){min(which(x>0))})
x$H3K27me3 <- apply(x[, 109:120], 1, function(x){min(which(x>0))})
x <- x[, 121:130]

# 10. convert time-point to real number of hours
hours2 <- c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168)


# 11. compute amount of anticipation / delay in hours
# between expression and each mark
x2 <- apply(x, 2, function(x){hours2[x]})
rownames(x2) <- rownames(x)
x2 <- as.data.frame(x2)
x2 <- apply(x2[, c(1:2, 4:10)], 2, function(x){x - x2$expression_tp})
x2[is.na(x2)] <- 170
x2 <- reshape2::melt(x2)
colnames(x2) <- c("gene_id", "mark", "value")

x2.promoters <- x2
# Proximal enhancers ----
rm(list=setdiff(ls(), c("expression.tp", "x2.promoters", "x.promoters")))


# 2. set working directory
setwd("../enhancers/proximal.enhancers.first.TSS/all.marks/")


# 3. read set of 257 activated genes
x <- read.table("../../../analysis/all.marks/expression/257.notExpressed.0h.txt",
                stringsAsFactors = F)
colnames(x) <- "gene_id"


# 4. the histone marks we analyze
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac","H3K4me3", "H3K36me3",
           "H4K20me1", "H3K9me3", "H3K27me3")


# 5. list of time-points
hours <- c("H000", "H003", "H006", "H009","H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")


# 6. for every mark, read dataframe with info 
# of presence/absence of peaks across time-points
# for the 257 activated genes
for ( i in 1:9) {
  
  # 6.1. read dataframe with 
  # 1st coln = gene_id
  # 2nd coln = time-points with peaks
  tmp <- read.table(paste0(marks[i], 
                           "/QN.merged/ant.del.analysis/", 
                           marks[i],
                           ".peaks.dynamics.tsv"), 
                    stringsAsFactors = F)
  
  colnames(tmp) <- c("gene_id", "tp")
  
  # 6.2. fix gene_id (aka remove everything after dot)
  tmp$gene_id <- gsub("\\..*","", tmp$gene_id)
  
  # 6.3. create empty dataframe
  m <- data.frame(stringsAsFactors = F)
  
  # 6.4. for each of the 257 activated genes
  for (gene in x$gene_id) {
    
    # 6.4.1. create a vector of 0s (n=12)
    v <- rep(0, 12)
    if (gene %in% tmp$gene_id) {
      
      # 6.4.2. retrieve the time-points in which the gene has a peak
      stn <- as.character(unlist(strsplit(tmp[tmp$gene_id == gene, "tp"], ";")))
      
      # 6.4.3. assign v=i only to the first time-point in which we have a peak
      v[min(which(hours %in% stn))] <- i
      
    }
    m <- rbind(m, v)
    
  }
  
  colnames(m) <- paste(hours, marks[i], sep="_")
  x <- cbind(x, m)
  
}


# 7. add rownames and remove gene_id column
rownames(x) <- x$gene_id
x$gene_id <- NULL


# 8. reorder columns in x according to order of marks
x <- x[, paste(rep(hours, 9), rep(c("H3K4me1", "H3K4me2",
                                     "H3K27ac", "H3K9ac", 
                                     "H3K4me3",
                                     "H3K36me3", "H4K20me1", 
                                     "H3K9me3", "H3K27me3"), each = 12), sep = "_")]


# 9. retrieve time-point of mark deposition / beginning of expression
x$H3K4me1 <- apply(x[, 1:12], 1, function(x){min(which(x>0))})
x$H3K4me2 <- apply(x[, 13:24], 1, function(x){min(which(x>0))})
x$H3K27ac <- apply(x[, 25:36], 1, function(x){min(which(x>0))})
x$H3K9ac <- apply(x[, 37:48], 1, function(x){min(which(x>0))})
x$H3K4me3 <- apply(x[, 49:60], 1, function(x){min(which(x>0))})
x$H3K36me3 <- apply(x[, 61:72], 1, function(x){min(which(x>0))})
x$H4K20me1 <- apply(x[, 73:84], 1, function(x){min(which(x>0))})
x$H3K9me3 <- apply(x[, 85:96], 1, function(x){min(which(x>0))})
x$H3K27me3 <- apply(x[, 97:108], 1, function(x){min(which(x>0))})
x <- x[, 109:117] 

stopifnot(identical(rownames(x), rownames(expression.tp)))
x <- merge(x, expression.tp, by = 0)
# we need to add +1, as in the code that we source 
# we were subtracting -1 to prepare for the heatmaps
x$expression_tp <- x$expression_tp + 1  
rownames(x) <- x$Row.names
x$Row.names <- NULL

# 10. convert time-point to real number of hours
hours2 <- c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168)


# 11. compute amount of anticipation / delay in hours
# between expression and each mark
x2 <- apply(x, 2, function(x){hours2[x]})
rownames(x2) <- rownames(x)
x2 <- as.data.frame(x2)
x2 <- apply(x2[, 1:9], 2, function(x){x - x2$expression_tp})
x2[is.na(x2)] <- 170
x2 <- melt(x2)
colnames(x2) <- c("gene_id", "mark", "value")

x2.proxEnhancers <- x2

# Distal enhancers ----
setwd("../../../analysis/enhancers/ENCODE-rE2G.blood/all.marks/")

# 3. read set of 257 activated genes
x <- read.table("../../../analysis/all.marks/expression/257.notExpressed.0h.txt",
                stringsAsFactors = F)
colnames(x) <- "gene_id"

# 4. the histone marks we analyze
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac","H3K4me3", "H3K36me3",
           "H4K20me1", "H3K9me3", "H3K27me3")

# 5. list of time-points
hours <- c("H000", "H003", "H006", "H009","H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")


# 6. for every mark, read dataframe with info 
# of presence/absence of peaks across time-points
# for the 257 activated genes

# 7.0. this list will store the dataframes of presence/absence of peaks
# for every histone mark
# rownames might be different among dfs because not all enhancers are marked by all marks
y <- list()

for ( i in 1:9) {
  
  # 7.1. read dataframe with 
  # 1st coln = gene_id
  # 2nd coln = time-points with peaks
  tmp <- read.table(paste0(marks[i], 
                           "/QN.merged/ant.del.analysis/", 
                           marks[i],
                           ".peaks.dynamics.tsv"), 
                    stringsAsFactors = F)
  
  tmp <- separate(tmp, "V1",  c("gene_id", "enhancer_id"),  sep = "chr")
  tmp$enhancer_id <- paste0("chr", tmp$enhancer_id)
  colnames(tmp)[3] <- c("tp")
  
  m <- data.frame(stringsAsFactors = F)
  
  # 7.2. keep from tmp only the 257 activated genes
  tmp <- tmp[tmp$gene_id %in% x$gene_id, ]
  
  # if there is at least one enhancer-gene pair marked by a given histone mark 
  if (nrow(tmp) > 0) { 
    
    # 7.3. add column of pair_id
    # this is needed when dfs for different marks are merged,
    # to ensure that the same gene-enhancer pair is correctly merged
    # across multiple marks
    tmp$pair_id <- paste0(tmp$gene_id, tmp$enhancer_id)
    
    for (k in 1:nrow(tmp)) {
      
      # 7.4. record the first time-point when the enhancer is marked, all other tps will be zeros 
      v <- rep(0, 12)
      stn <- as.character(unlist(strsplit(tmp[k, "tp"], ";")))
      v[min(which(hours %in% stn))] <- i
      m <- rbind(m, v)
      
    }
    
    colnames(m) <- paste(hours, marks[i], sep="_")
    rownames(m) <- tmp$pair_id
    
    # 7.5. otherwise, set m as a 1-row matrix of zeros  
    # we need m to be non-empty, since later we loop across all elements of
    # the y list to merge dfs of different histone marks
  } else {
    
    m <- matrix(rep(0, 12), byrow = T, nrow = 1)
    colnames(m) <- paste(hours, marks[i], sep="_")
    m <- as.data.frame(m)
  }
  
  y[[i]] <- m
  
  
}


# 8. collect all possible pairs of gene_enhancers marked by any mark
# in the case of ABC.genewise.blood.closest, we can't have more than one enhancer per gene
all.pairs <- c()

for ( i in 1:length(y)) {
  
  all.pairs <- c(all.pairs, rownames(y[[i]]))
  
}

all.pairs <- unique(all.pairs)


# 9. create a unique dataframe across all marks
x.enhancers <- data.frame(pair_id = all.pairs, stringsAsFactors = F)
for ( i in 1:length(y)) {
  
  y[[i]]$pair_id <- rownames(y[[i]])
  x.enhancers <- merge(x.enhancers, y[[i]], all.x = T, by = "pair_id")
  
}


# 10. reorder x.enhancers columns
x.enhancers <- x.enhancers[, c("pair_id", paste(rep(hours, 9), rep(c("H3K4me1", 
                                                                      "H3K4me2", "H3K27ac", 
                                                                      "H3K9ac", "H3K4me3",
                                                                      "H3K36me3", "H4K20me1", 
                                                                      "H3K9me3", "H3K27me3"), each = 12), sep = "_"))]


# 11. substitute NAs with 0s
x.enhancers[is.na(x.enhancers)] <- 0
x.enhancers <- x.enhancers[grep("ENSG", x.enhancers$pair_id), ]


# 12. retrieve gene_ids of genes with an enhancer
x.enhancers <- separate(x.enhancers, "pair_id",  c("gene_id", "enhancer_id"),  sep = "chr")
x.enhancers$enhancer_id <- NULL


# 13. get list of genes with one enhancer from Vasilis' list
# genes.NA contains genes that don't have enhancers, they will be set as NAs in the final x matrix
genes.pairs <- read.delim("../interacting_pairs_filtered_100bp.merged_closest.bed",
                          h=F)
genes.pairs <- genes.pairs$V4
genes.pairs <- gsub("\\c.*","", genes.pairs)
genes.NA <- setdiff(x$gene_id, genes.pairs) 


# 14. add to x:
# - the genes unmarked at enhancers for all marks (they will be all zeros)
# - the genes that don't have enhancers (they will be all NAs)
x <- merge(x, x.enhancers, by = "gene_id", all = T)
x[is.na(x)] <- 0
rownames(x) <- x$gene_id
x$gene_id <- NULL
x[rownames(x) %in% genes.NA, ] <- NA


# 15. retrieve time-point of mark deposition / beginning of expression
x$H3K4me1 <- apply(x[, 1:12], 1, function(x){min(which(x>0))})
x$H3K4me2 <- apply(x[, 13:24], 1, function(x){min(which(x>0))})
x$H3K27ac <- apply(x[, 25:36], 1, function(x){min(which(x>0))})
x$H3K9ac <- apply(x[, 37:48], 1, function(x){min(which(x>0))})
x$H3K4me3 <- apply(x[, 49:60], 1, function(x){min(which(x>0))})
x$H3K36me3 <- apply(x[, 61:72], 1, function(x){min(which(x>0))})
x$H4K20me1 <- apply(x[, 73:84], 1, function(x){min(which(x>0))})
x$H3K9me3 <- apply(x[, 85:96], 1, function(x){min(which(x>0))})
x$H3K27me3 <- apply(x[, 97:108], 1, function(x){min(which(x>0))})
x <- x[, 109:117] 

stopifnot(identical(rownames(x), rownames(expression.tp)))
x <- merge(x, expression.tp, by = 0)
# we need to add +1, as in the code that we source 
# we were subtracting -1 to prepare for the heatmaps
x$expression_tp <- x$expression_tp + 1  
rownames(x) <- x$Row.names
x$Row.names <- NULL

# 16. convert time-point to real number of hours
hours2 <- c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168)


# 17. compute amount of anticipation / delay in hours
# between expression and each mark
x2 <- apply(x, 2, function(x){hours2[x]})
rownames(x2) <- rownames(x)
x2 <- as.data.frame(x2)
x2 <- apply(x2[, 1:9], 2, function(x){x - x2$expression_tp})
x2[is.na(x2)] <- 170
x2 <- melt(x2)
colnames(x2) <- c("gene_id", "mark", "value")

x2.distEnhancers <- x2

# Integrate data from promoters, proximal and distal enhancers ----
rm(list=setdiff(ls(), c("x2.promoters", "x2.proxEnhancers", "x2.distEnhancers")))

x2.promoters$region <- "promoters"
x2.proxEnhancers$region <- "proximal enhancers"
x2.distEnhancers$region <- "distal enhancers"

# 2. merge the three types of regions in a unique df
x2 <- rbind(x2.promoters, x2.proxEnhancers, x2.distEnhancers)
x2$combo <- paste(x2$mark, x2$region, sep = "_")

x2 <- x2[complete.cases(x2), ]
x2$region <- factor(x2$region, levels = c("promoters", "proximal enhancers", "distal enhancers"))
x2 <- x2[x2$mark %in% c("H3K4me1",
                        "H3K4me2",
                        "H3K27ac",
                        "H3K9ac",
                        "H3K4me3",
                        "H3K36me3",
                        "H4K20me1",
                        "H3K9me3",
                        "H3K27me3"), ]
x2$mark <- factor(x2$mark,
                  levels = rev(c("H3K4me1", "H3K4me2",
                                 "H3K27ac", "H3K9ac",
                                 "H3K4me3",
                                 "H3K36me3",
                                 "H4K20me1",
                                 "H3K9me3",
                                 "H3K27me3")))

# 5. define color palette
palette <- c("H3K9ac" = "#af4d85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#d199b9",
             "H3K27me3" = "#1d2976",
             "H3K9me3" = "#a7add4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#4c7027",
             "H3K4me1" = "#e5ab00",
             "H3K4me2" = "#a67c00")


# 6. make plot
pdf("fig_4g.pdf", 
    height = 4, width = 7, useDingbats = F)
ggplot(data = x2) +
  geom_violin(aes(x=mark, y=value, fill=mark),
              alpha=.5, width = 1, color = "white") +
  stat_summary(mapping = aes(x = mark, y = value),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, size = .2) +
  # geom_boxplot(data = x2.unmarked.H3K27me3, aes(x=mark, y=value, fill=mark),
  #              width=0.2, alpha = .8, outlier.shape = NA) +
  # geom_text(data = text.df, aes(y=pos, x=Var1, label = paste0(Freq, "%"))) +
  facet_grid(.~region) +
  guides(fill='none') +
  theme_bw() +
  coord_flip() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=12, angle = 45, vjust = .95, hjust = .95),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        panel.border = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        plot.title = element_text(size = 12, hjust = .5)) +
  scale_fill_manual(values = palette) +
  geom_hline(yintercept = 0, color = "#525252", linetype="dashed") +
  ylab("time of histone marking") +
  scale_y_continuous(labels = c("100 h\nbefore", "gene\nactivation", "100 h\nafter", "never\nmarked"),
                     breaks = c(-100, 0, 100, 170))
dev.off()
