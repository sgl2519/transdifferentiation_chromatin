#************
# LIBRARIES *
#************

library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(RColorBrewer)

#********
# BEGIN *
#********

# 1. set wd
setwd("../analysis/all.marks")

# 2. read list of 257 not expressed genes
x <- read.table("expression/257.notExpressed.0h.txt",
                stringsAsFactors = F)
colnames(x) <- "gene_id"

# 3. marks 
marks <- c("H3K27ac", "H3K9ac","H3K4me3", "H3K4me1", "H3K4me2", "H3K36me3",
           "H4K20me1", "H3K9me3", "H3K27me3", "cEBPa")

# 4. time-points
hours <- c("H000", "H003", "H006", "H009","H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")
hours2 <- c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168)
names(hours2) <- hours

# 5. retrieve time-point of peak appearance of 
# for each histone mark

for ( i in c(1:10)) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/ant.del.analysis/", marks[i],
                             ".peaks.dynamics.tsv"), stringsAsFactors = F)
      
  colnames(tmp) <- c("gene_id", "tp")
  tmp$gene_id <- gsub("\\..*","", tmp$gene_id)
  tmp$tp <- gsub("\\;.*","", tmp$tp)
  colnames(tmp)[2] <- marks[i]
  x <- merge(x, tmp, by = "gene_id", all.x = T)
  
}

# 6. add time-point for expression (aka the time-point in which exp is > 1 TPM)
expression.matrix <- read.table("H3K4me3/QN.merged/expression.matrix.tsv", h=T, sep="\t")
expression.matrix <- expression.matrix[rownames(expression.matrix) %in% x$gene_id, ]

m <- data.frame(stringsAsFactors = F)

for (gene in x$gene_id) {
  
  v <- rep(0, 12)
  stn <- as.numeric(expression.matrix[gene, ])
  
  v[which(stn > 1)] <- 1
  
  m <- rbind(m, v)
}

x$expression <- hours[apply(m, 1, function(x){min(which(x>0))})]


# 7. remove genes marked by K27me3 at 0 hours
x <- x[(x$H3K27me3 == "H000" | x$H3K27me3 == "H003") & !(is.na(x$H3K27me3)) & !(is.na(x$cEBP)), ]

# 8. change the hours format to integers
for ( i in 2:ncol(x) ){
  
  x[, i] <- hours2[x[, i]]
  
}

# 9. compute pair-wise comparisons between marks/expression
# aka matrix for marks on rows anticipating marks on columns
df.anticipate <- data.frame(stringsAsFactors = F)
for (i in colnames(x)[2:ncol(x)]) {
  
  for ( j in colnames(x)[2:ncol(x)]) {
    
    if (i != j) {
      
      tmp <- as.data.frame(x[, c(i, j)])
      tmp <- tmp[!is.na(tmp[, 1]), ]
      tmp[is.na(tmp)] <- 170
      tmp$diff <- tmp[, 1] - tmp[, 2]
      
      tmp2 <- data.frame(n = nrow(tmp[tmp$diff < 0, ]),
                         mark1 = i,
                         mark2 = j)
      
      tot <- nrow(x) - nrow(x[is.na(x[, i]) & is.na(x[, j]), ])
      tmp2$n <- round(tmp2$n/tot*100, 1)
      
      df.anticipate <- rbind(df.anticipate, tmp2)
      
    }
    
  }
  
}

df.anticipate <- reshape2::dcast(df.anticipate, formula=mark1~mark2, value.var = "n")
rownames(df.anticipate) <- df.anticipate$mark1
df.anticipate$mark1 <- NULL
sorted.marks <- c("H3K4me1", "H3K4me2", "H3K4me3", "cEBPa",
                  "H4K20me1", "H3K9ac", "expression",
                  "H3K27ac", "H3K36me3")
df.anticipate <- df.anticipate[sorted.marks, sorted.marks]
rownames(df.anticipate) <- colnames(df.anticipate)

# 10. compute pair-wise comparisons between marks/expression
# aka matrix for marks on rows co-occurring with marks on columns
df.concomitant <- data.frame(stringsAsFactors = F)
for (i in colnames(x)[2:ncol(x)]) {
  
  for ( j in colnames(x)[2:ncol(x)]) {
    
    if (i != j) {
      
      tmp <- as.data.frame(x[, c(i, j)])
      tmp <- tmp[!is.na(tmp[, 1]), ]
      tmp[is.na(tmp)] <- 170
      tmp$diff <- tmp[, 1] - tmp[, 2]
      
      tmp2 <- data.frame(n = nrow(tmp[tmp$diff == 0, ]),
                         mark1 = i,
                         mark2 = j)
      
      tot <- nrow(x) - nrow(x[is.na(x[, i]) & is.na(x[, j]), ])
      tmp2$n <- round(tmp2$n/tot*100, 1)
      
      df.concomitant <- rbind(df.concomitant, tmp2)
      
      
    }
    
  }
  
}

df.concomitant <- reshape2::dcast(df.concomitant, formula=mark1~mark2, value.var = "n")
rownames(df.concomitant) <- df.concomitant$mark1
df.concomitant$mark1 <- NULL
df.concomitant <- df.concomitant[sorted.marks, sorted.marks]
rownames(df.concomitant) <- colnames(df.concomitant)

# Find the best order
comb <- combinat::permn(c(colnames(df.anticipate)))
s <- c()
for ( c in comb ) {
  tmp <- df.anticipate[unlist(c), unlist(c)]
  s <- c(s, 
         (sum(tmp[upper.tri(tmp)]) - sum(tmp[lower.tri(tmp)]))
  )
}

pdf("fig_S10d_right.matrix.pdf",
    width = 5, height = 5)
pheatmap(as.matrix(df.anticipate[unlist(comb[which.max(s)]), unlist(comb[which.max(s)])]), 
         cluster_rows = F, cluster_cols = F, display_numbers = T,
         border_color = NA, fontsize = 15, legend = F, 
         number_format = "%.1f",
         color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(19),
         breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90))
dev.off()
