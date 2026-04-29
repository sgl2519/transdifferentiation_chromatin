#************
# LIBRARIES *
#************

library(ggplot2)
library(pheatmap)


#************
# FUNCTIONS *
#************

f.anticipate <- function(mark) {
  
  
  x.sub <- x[, c(paste0(mark, "_promoter"), paste0(mark, "_pEnh"), paste0(mark, "_dEnh"))]
  df.anticipate <- data.frame(stringsAsFactors = F)
  for (i in colnames(x.sub)[1:ncol(x.sub)]) {
    
    for ( j in colnames(x.sub)[1:ncol(x.sub)]) {
      
      if (i != j) {
        
        tmp <- as.data.frame(x.sub[, c(i, j)])
        tmp <- tmp[!is.na(tmp[, 1]), ]
        tmp[is.na(tmp)] <- 170
        tmp$diff <- tmp[, 1] - tmp[, 2]
        
        tmp2 <- data.frame(n = nrow(tmp[tmp$diff < 0, ]),
                           mark1 = i,
                           mark2 = j)
        
        tot <- nrow(x) - nrow(x[is.na(x.sub[, i]) & is.na(x.sub[, j]), ])
        tmp2$n <- round(tmp2$n/tot*100, 1)
        
        df.anticipate <- rbind(df.anticipate, tmp2)
        
      }
      
    }
    
  }
  
  df.anticipate <- reshape2::dcast(df.anticipate, formula=mark1~mark2, value.var = "n")
  rownames(df.anticipate) <- df.anticipate$mark1
  df.anticipate$mark1 <- NULL
  sorted.marks <- c(c(paste0(marks[k], "_promoter"), paste0(marks[k], "_pEnh"), paste0(marks[k], "_dEnh")))
  df.anticipate <- df.anticipate[sorted.marks, sorted.marks]
  colnames(df.anticipate) <- c("prom.", "pEnh", "dEnh")
  rownames(df.anticipate) <- colnames(df.anticipate)
  
  p <- pheatmap::pheatmap(as.matrix(df.anticipate), 
                cluster_rows = F, cluster_cols = F, display_numbers = T,
                border_color = "black", fontsize = 15, legend = F,
                na_col = "white",
                cellwidth = 50, cellheight = 50,
                number_format = "%.1f",
                main = mark,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(15),
                breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70))
  
  return(p$gtable)
  
}
f.concomitant <- function(mark) {
  
  x.sub <- x[, c(paste0(mark, "_promoter"), paste0(mark, "_pEnh"), paste0(mark, "_dEnh"))]
  df.concomitant <- data.frame(stringsAsFactors = F)
  for (i in colnames(x.sub)[1:ncol(x.sub)]) {
    
    for ( j in colnames(x.sub)[1:ncol(x.sub)]) {
      
      if (i != j) {
        
        tmp <- as.data.frame(x.sub[, c(i, j)])
        tmp <- tmp[!is.na(tmp[, 1]), ]
        tmp[is.na(tmp)] <- 170
        tmp$diff <- tmp[, 1] - tmp[, 2]
        
        tmp2 <- data.frame(n = nrow(tmp[tmp$diff == 0, ]),
                           mark1 = i,
                           mark2 = j)
        
        tot <- nrow(x) - nrow(x[is.na(x.sub[, i]) & is.na(x.sub[, j]), ])
        tmp2$n <- round(tmp2$n/tot*100, 1)
        
        df.concomitant <- rbind(df.concomitant, tmp2)
        
      }
      
    }
    
  }
  
  df.concomitant <- reshape2::dcast(df.concomitant, formula=mark1~mark2, value.var = "n")
  rownames(df.concomitant) <- df.concomitant$mark1
  df.concomitant$mark1 <- NULL
  sorted.marks <- c(c(paste0(marks[k], "_promoter"), paste0(marks[k], "_pEnh"), paste0(marks[k], "_dEnh")))
  df.concomitant <- df.concomitant[sorted.marks, sorted.marks]
  colnames(df.concomitant) <- c("prom.", "pEnh", "dEnh")
  rownames(df.concomitant) <- colnames(df.concomitant)
  
  p <- pheatmap::pheatmap(as.matrix(df.concomitant), 
                cluster_rows = F, cluster_cols = F, display_numbers = T,
                border_color = "black", fontsize = 15, legend = F,
                na_col = "white",
                cellwidth = 50, cellheight = 50,
                number_format = "%.1f",
                main = mark,
                color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(15),
                breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70))
  
  return(p$gtable)
  
}


#********
# BEGIN *
#********


# Generate matrix of peak appearance for promoters ----
setwd("../analysis/all.marks")

# 2. read list of 257 not expressed genes
x <- read.table("expression/257.notExpressed.0h.txt",
                stringsAsFactors = F)
colnames(x) <- "gene_id"

# 3. marks 
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3")

# 4. time-points
hours <- c("H000", "H003", "H006", "H009","H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")
hours2 <- c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168)
names(hours2) <- hours

# 5. retrieve time-point of peak appearance of 
# for each histone mark

for ( i in c(1:5)) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/ant.del.analysis/", marks[i],
                           ".peaks.dynamics.tsv"), stringsAsFactors = F)
  colnames(tmp) <- c("gene_id", "tp")
  tmp$gene_id <- gsub("\\..*","", tmp$gene_id)
  tmp$tp <- gsub("\\;.*","", tmp$tp)
  colnames(tmp)[2] <- marks[i]
  x <- merge(x, tmp, by = "gene_id", all.x = T)
  
}

# 6. add time-point for expression (aka the time-point in which exp is > 1 TPM)
expression.matrix <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
expression.matrix <- expression.matrix[rownames(expression.matrix) %in% x$gene_id, ]

m <- data.frame(stringsAsFactors = F)

for (gene in x$gene_id) {
  
  v <- rep(0, 12)
  stn <- as.numeric(expression.matrix[gene, ])
  
  v[which(stn > 1)] <- 1
  
  m <- rbind(m, v)
}

x$expression <- hours[apply(m, 1, function(x){min(which(x>0))})]


# 8. change the hours format to integers
for ( i in 2:ncol(x) ){
  
  x[, i] <- hours2[x[, i]]
  
}

x.promoter <- x
rownames(x.promoter) <- x.promoter$gene_id
x.promoter$gene_id <- NULL
colnames(x.promoter) <- paste0(colnames(x.promoter), '_promoter')

# Generate matrix of peak appearance for pEnh ----
setwd("../enhancers/proximal.enhancers.first.TSS/all.marks/")

# 2. read list of 257 not expressed genes
x <- read.table("../../../all.marks/expression/257.notExpressed.0h.txt",
                stringsAsFactors = F)
colnames(x) <- "gene_id"

# 3. marks 
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3")

# 4. time-points
hours <- c("H000", "H003", "H006", "H009","H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")
hours2 <- c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168)
names(hours2) <- hours

# 5. retrieve time-point of peak appearance of 
# for each histone mark

for ( i in c(1:5)) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/ant.del.analysis/", marks[i],
                           ".peaks.dynamics.tsv"), stringsAsFactors = F)
  colnames(tmp) <- c("gene_id", "tp")
  tmp$gene_id <- gsub("\\..*","", tmp$gene_id)
  tmp$tp <- gsub("\\;.*","", tmp$tp)
  colnames(tmp)[2] <- marks[i]
  x <- merge(x, tmp, by = "gene_id", all.x = T)
  
}

# 6. add time-point for expression (aka the time-point in which exp is > 1 TPM)
expression.matrix <- read.table("../../../all.marks/expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
expression.matrix <- expression.matrix[rownames(expression.matrix) %in% x$gene_id, ]

m <- data.frame(stringsAsFactors = F)

for (gene in x$gene_id) {
  
  v <- rep(0, 12)
  stn <- as.numeric(expression.matrix[gene, ])
  
  v[which(stn > 1)] <- 1
  
  m <- rbind(m, v)
}

x$expression <- hours[apply(m, 1, function(x){min(which(x>0))})]


# 8. change the hours format to integers
for ( i in 2:ncol(x) ){
  
  x[, i] <- hours2[x[, i]]
  
}

x.pEnh <- x
rownames(x.pEnh) <- x.pEnh$gene_id
x.pEnh$gene_id <- NULL
colnames(x.pEnh) <- paste0(colnames(x.pEnh), '_pEnh')

# Generate matrix of peak appearance for dEnh ----
setwd("../../../enhancers/ENCODE-rE2G.blood/all.marks/")

# 2. read list of 257 not expressed genes
x <- read.table("../../../all.marks/expression/257.notExpressed.0h.txt",
                stringsAsFactors = F)
colnames(x) <- "gene_id"

# 3. marks 
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3")

# 4. time-points
hours <- c("H000", "H003", "H006", "H009","H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")
hours2 <- c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168)
names(hours2) <- hours

# 5. retrieve time-point of peak appearance of 
# for each histone mark

for ( i in c(1:5)) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/ant.del.analysis/", marks[i],
                           ".peaks.dynamics.tsv"), stringsAsFactors = F)
  tmp$V1 <- gsub('c.*', '', tmp$V1)
  colnames(tmp) <- c("gene_id", "tp")
  tmp$gene_id <- gsub("\\..*","", tmp$gene_id)
  tmp$tp <- gsub("\\;.*","", tmp$tp)
  colnames(tmp)[2] <- marks[i]
  x <- merge(x, tmp, by = "gene_id", all.x = T)
  
}

# 6. add time-point for expression (aka the time-point in which exp is > 1 TPM)
expression.matrix <- read.table("../../../all.marks/expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
expression.matrix <- expression.matrix[rownames(expression.matrix) %in% x$gene_id, ]

m <- data.frame(stringsAsFactors = F)

for (gene in x$gene_id) {
  
  v <- rep(0, 12)
  stn <- as.numeric(expression.matrix[gene, ])
  
  v[which(stn > 1)] <- 1
  
  m <- rbind(m, v)
}

x$expression <- hours[apply(m, 1, function(x){min(which(x>0))})]


# 8. change the hours format to integers
for ( i in 2:ncol(x) ){
  
  x[, i] <- hours2[x[, i]]
  
}

x.dEnh <- x
rownames(x.dEnh) <- x.dEnh$gene_id
x.dEnh$gene_id <- NULL
colnames(x.dEnh) <- paste0(colnames(x.dEnh), '_dEnh')

# Plot matrices
x <- cbind(x.promoter, x.pEnh, x.dEnh)

# 3. compute differences
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3")

lop <- list()
for(k in 1:5) {
  
  lop[[k]] <- f.anticipate(mark = marks[k])
}
dev.off()

for(k in 1:5) {
  
  lop[[k+5]] <- f.concomitant(mark = marks[k])
}
dev.off()

pdf("fig_S11b.pdf",
    width = 15, height = 7)
plot_grid(plotlist = lop, nrow=2, ncol=5)
dev.off()
