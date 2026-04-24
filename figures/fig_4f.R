library(seriation)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)

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
           "H4K20me1", "H3K9me3", "H3K27me3", "cEBP")

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


# 9. the seven active histone marks
marks <- marks[-c(10)]

# 10. columns corresponding to the seven active marks
sorted.tps <- c(1, # H3K4me1
                13, # H3K4me2
                37, # H3K27ac
                49, # H3K9ac
                61, # H3K4me3
                73, # H3K36me3
                85, # H4K20me1
                97, # H3K9me3
                109)# H3K27me3

# 11. select genes unmarked by H3K27me3 at 0-3h
unmarked.H3K27me3 <- rownames(x[x$H000_H3K27me3 == 0 & x$H003_H3K27me3 == 0, ])
marked.H3K27me3 <- rownames(x[x$H000_H3K27me3 == 9 | x$H003_H3K27me3 == 9, ])


# 12. for each gene, retrieve the first time-point of expression > 1 TPM 
x$expression_tp <- apply(x[, 25:36], 1, function(x){min(which(x > 0))})
# 12.1. subtract -1 to expression_tp to make easier comparison when preparing the matrices below
x$expression_tp <- x$expression_tp - 1 
# 12.2. save dataframe to a new variable - this is needed when working on enhancers
expression.tp <- x[, "expression_tp", drop = F]


#----------------------------------------------------
# 13. subset genes based on time-point of expression
#----------------------------------------------------

# 13.1. genes expressed at 3 hours
x.3h <- x[x$expression_tp == 1, ]
x.3h <- cbind(x.3h[, sorted.tps], # 0h 
              x.3h[, sorted.tps+1], # 3h
              (x.3h[, sorted.tps+2] + x.3h[, sorted.tps+3] + x.3h[, sorted.tps+4] + 
                 x.3h[, sorted.tps+5] + x.3h[, sorted.tps+6] + x.3h[, sorted.tps+7] + 
                 x.3h[, sorted.tps+8] + x.3h[, sorted.tps+9] + x.3h[, sorted.tps+10] + x.3h[, sorted.tps+11])) # 6-168h
colnames(x.3h) <- c(paste0("before_", marks),
                    paste0("expression_", marks),
                    paste0("after_", marks))

# 13.2. genes expressed at 6 hours
x.6h <- x[x$expression_tp == 2, ]
x.6h <- cbind((x.6h[, sorted.tps] + x.6h[, sorted.tps+1]), # 0-3h 
              x.6h[, sorted.tps+2], # 6h
              (x.6h[, sorted.tps+3] + x.6h[, sorted.tps+4] + x.6h[, sorted.tps+5] + 
                 x.6h[, sorted.tps+6] + x.6h[, sorted.tps+7] + x.6h[, sorted.tps+8] +
                 x.6h[, sorted.tps+9] + x.6h[, sorted.tps+10] + x.6h[, sorted.tps+11])) # 9-168h
colnames(x.6h) <- c(paste0("before_", marks),
                    paste0("expression_", marks),
                    paste0("after_", marks))

# 13.3. genes expressed at 9 hours
x.9h <- x[x$expression_tp == 3, ]
x.9h <- cbind((x.9h[, sorted.tps] + x.9h[, sorted.tps+1] + x.9h[, sorted.tps+2]), #0-6h
              x.9h[, sorted.tps+3], #9h
              (x.9h[, sorted.tps+4] + x.9h[, sorted.tps+5] + x.9h[, sorted.tps+6] + 
                 x.9h[, sorted.tps+7] + x.9h[, sorted.tps+8] + x.9h[, sorted.tps+9] + 
                 x.9h[, sorted.tps+10] + x.9h[, sorted.tps+11])) # 12-168h
colnames(x.9h) <- c(paste0("before_", marks),
                    paste0("expression_", marks),
                    paste0("after_", marks))

# 13.4. genes expressed at 12 hours
x.12h <- x[x$expression_tp == 4, ]
x.12h <- cbind((x.12h[, sorted.tps] + x.12h[, sorted.tps+1] + x.12h[, sorted.tps+2] + x.12h[, sorted.tps+3]), #0-9h
               x.12h[, sorted.tps+4], #12h
               (x.12h[, sorted.tps+5] + x.12h[, sorted.tps+6] + x.12h[, sorted.tps+7] + 
                  x.12h[, sorted.tps+8] + x.12h[, sorted.tps+9] + x.12h[, sorted.tps+10] + x.12h[, sorted.tps+11])) #18-168h
colnames(x.12h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 13.5. genes expressed at 18 hours
x.18h <- x[x$expression_tp == 5, ]
x.18h <- cbind((x.18h[, sorted.tps] + x.18h[, sorted.tps+1] + x.18h[, sorted.tps+2] + 
                  x.18h[, sorted.tps+3] + x.18h[, sorted.tps+4]), #0-12h
               x.18h[, sorted.tps+5], #18h
               (x.18h[, sorted.tps+6] + x.18h[, sorted.tps+7] + x.18h[, sorted.tps+8] +
                  x.18h[, sorted.tps+9] + x.18h[, sorted.tps+10] + x.18h[, sorted.tps+11])) #24-168h
colnames(x.18h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 13.6. genes expressed at 24 hours
x.24h <- x[x$expression_tp == 6, ]
x.24h <- cbind((x.24h[, sorted.tps] + x.24h[, sorted.tps+1] + x.24h[, sorted.tps+2] + 
                  x.24h[, sorted.tps+3] + x.24h[, sorted.tps+4] + x.24h[, sorted.tps+5]), #0-18h
               x.24h[, sorted.tps+6], #24h
               (x.24h[, sorted.tps+7] + x.24h[, sorted.tps+8] + x.24h[, sorted.tps+9] + 
                  x.24h[, sorted.tps+10] + x.24h[, sorted.tps+11])) #36-168h
colnames(x.24h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 13.7. genes expressed at 36 hours
x.36h <- x[x$expression_tp == 7, ]
x.36h <- cbind((x.36h[, sorted.tps] + x.36h[, sorted.tps+1] + x.36h[, sorted.tps+2] + 
                  x.36h[, sorted.tps+3] + x.36h[, sorted.tps+4] + x.36h[, sorted.tps+5] +
                  x.36h[, sorted.tps+6]), #0-24h
               x.36h[, sorted.tps+7], #36h
               (x.36h[, sorted.tps+8] + x.36h[, sorted.tps+9] + x.36h[, sorted.tps+10] + x.36h[, sorted.tps+11])) #48-168h
colnames(x.36h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 13.8. genes expressed at 48 hours
x.48h <- x[x$expression_tp == 8, ]
x.48h <- cbind((x.48h[, sorted.tps] + x.48h[, sorted.tps+1] + x.48h[, sorted.tps+2] + 
                  x.48h[, sorted.tps+3] + x.48h[, sorted.tps+4] + x.48h[, sorted.tps+5] + 
                  x.48h[, sorted.tps+6] + x.48h[, sorted.tps+7]), #0-36h
               x.48h[, sorted.tps+8], #48h
               (x.48h[, sorted.tps+9] + x.48h[, sorted.tps+10] + x.48h[, sorted.tps+11])) #72-168h
colnames(x.48h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 13.9. genes expressed at 72 hours
x.72h <- x[x$expression_tp == 9, ]
x.72h <- cbind((x.72h[, sorted.tps] + x.72h[, sorted.tps+1] + x.72h[, sorted.tps+2] + 
                  x.72h[, sorted.tps+3] + x.72h[, sorted.tps+4] + x.72h[, sorted.tps+5] + 
                  x.72h[, sorted.tps+6] + x.72h[, sorted.tps+7] + x.72h[, sorted.tps+8]), #0-48h
               x.72h[, sorted.tps+9], #72h
               (x.72h[, sorted.tps+10] + x.72h[, sorted.tps+11])) #120-168h 
colnames(x.72h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 13.10. genes expressed at 120 hours
x.120h <- x[x$expression_tp == 10, ]
x.120h <- cbind((x.120h[, sorted.tps] + x.120h[, sorted.tps+1] + x.120h[, sorted.tps+2] + 
                   x.120h[, sorted.tps+3] + x.120h[, sorted.tps+4] + x.120h[, sorted.tps+5] + 
                   x.120h[, sorted.tps+6] + x.120h[, sorted.tps+7] + x.120h[, sorted.tps+8] + x.120h[, sorted.tps+9]), #0-72h
                x.120h[, sorted.tps+10], #120h
                x.120h[, sorted.tps+11]) #168h 
colnames(x.120h) <- c(paste0("before_", marks),
                      paste0("expression_", marks),
                      paste0("after_", marks))

# 14. rbind all matrices into a unique dataframe
x <- rbind(x.3h, x.6h, x.9h, x.12h, x.18h, x.24h, x.36h, x.48h, x.72h, x.120h)

# 15. if a gene has already acquired a mark in a previous time-point,
# then it will have the mark also in the next time-points

x$expression_H3K4me1 <- ifelse(x$before_H3K4me1 > 0, 1, x$expression_H3K4me1)
x$after_H3K4me1 <- ifelse(x$expression_H3K4me1 > 0, 1, x$after_H3K4me1)

x$expression_H3K4me2 <- ifelse(x$before_H3K4me2 > 0, 2, x$expression_H3K4me2)
x$after_H3K4me2 <- ifelse(x$expression_H3K4me2 > 0, 2, x$after_H3K4me2)

x$expression_H3K27ac <- ifelse(x$before_H3K27ac > 0, 3, x$expression_H3K27ac)
x$after_H3K27ac <- ifelse(x$expression_H3K27ac > 0, 3, x$after_H3K27ac)

x$expression_H3K9ac <- ifelse(x$before_H3K9ac > 0, 4, x$expression_H3K9ac)
x$after_H3K9ac <- ifelse(x$expression_H3K9ac > 0, 4, x$after_H3K9ac)

x$expression_H3K4me3 <- ifelse(x$before_H3K4me3 > 0, 5, x$expression_H3K4me3)
x$after_H3K4me3 <- ifelse(x$expression_H3K4me3 > 0, 5, x$after_H3K4me3)

x$expression_H3K36me3 <- ifelse(x$before_H3K36me3 > 0, 6, x$expression_H3K36me3)
x$after_H3K36me3 <- ifelse(x$expression_H3K36me3 > 0, 6, x$after_H3K36me3)

x$expression_H4K20me1 <- ifelse(x$before_H4K20me1 > 0, 7, x$expression_H4K20me1)
x$after_H4K20me1 <- ifelse(x$expression_H4K20me1 > 0, 7, x$after_H4K20me1)

x$expression_H3K9me3 <- ifelse(x$before_H3K9me3 > 0, 8, x$expression_H3K9me3)
x$after_H3K9me3 <- ifelse(x$expression_H3K9me3 > 0, 8, x$after_H3K9me3)

x$expression_H3K27me3 <- ifelse(x$before_H3K27me3 > 0, 9, x$expression_H3K27me3)
x$after_H3K27me3 <- ifelse(x$expression_H3K27me3 > 0, 9, x$after_H3K27me3)


# 16. seriate matrix of genes unmarked by H3K27me3
x1 <- x[rownames(x) %in% unmarked.H3K27me3, ]
set.seed(1234)
# ser.obj1 <- seriate(as.matrix(x1[, c(1, 3, 10, 12, 19, 21)]), method = "PCA")
ser.obj1 <- seriate(as.matrix(x1[, c(1:5, 10:14, 19:23)]), method = "PCA")
x1.ser <- x1[get_order(ser.obj1, 1), ]


# 17. seriate matrix of genes marked by H3K27me3
x2 <- x[rownames(x) %in% marked.H3K27me3, ]
set.seed(1234)
# ser.obj2 <- seriate(as.matrix(x2[, c(1, 3, 10, 12, 19, 21)]), method = "PCA")
ser.obj2 <- seriate(as.matrix(x2[, c(1:5, 10:14, 19:23)]), method = "PCA")
x2.ser <- x2[get_order(ser.obj2, 1), ]

x.promoters <- x
x2.promoters <- x2

# Proximal enhancers ----
x1.ser.promoters <- x1.ser
x2.ser.promoters <- x2.ser
rm(list=setdiff(ls(), c("x1.ser.promoters", "x2.ser.promoters", "x2.promoters", "x.promoters", "expression.tp")))


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

# 9. ensure rownames are the same in x and expression.tp
stopifnot(identical(rownames(x), rownames(expression.tp)))

# 10. add to x the info on gene-expression time-point
x$expression_tp <- expression.tp$expression_tp


# 11. columns corresponding to the seven active marks
sorted.tps <- c(1, # H3K4me1
                13, # H3K4me2
                25, # H3K27ac
                37, # H3K9ac
                49, # H3K4me3
                61, # H3K36me3
                73, # H4K20me1
                85, # H3K9me3
                97) # H3K27me3

#marks <- marks[-9]


#----------------------------------------------------
# 12. subset genes based on time-point of expression
#----------------------------------------------------

# 12.1. genes expressed at 3 hours
x.3h <- x[x$expression_tp == 1, ]
x.3h <- cbind(x.3h[, sorted.tps], # 0h 
              x.3h[, sorted.tps+1], # 3h
              (x.3h[, sorted.tps+2] + x.3h[, sorted.tps+3] + x.3h[, sorted.tps+4] + 
                 x.3h[, sorted.tps+5] + x.3h[, sorted.tps+6] + x.3h[, sorted.tps+7] + 
                 x.3h[, sorted.tps+8] + x.3h[, sorted.tps+9] + x.3h[, sorted.tps+10] + x.3h[, sorted.tps+11])) # 6-168h
colnames(x.3h) <- c(paste0("before_", marks),
                    paste0("expression_", marks),
                    paste0("after_", marks))

# 12.2. genes expressed at 6 hours
x.6h <- x[x$expression_tp == 2, ]
x.6h <- cbind((x.6h[, sorted.tps] + x.6h[, sorted.tps+1]), # 0-3h 
              x.6h[, sorted.tps+2], # 6h
              (x.6h[, sorted.tps+3] + x.6h[, sorted.tps+4] + x.6h[, sorted.tps+5] + 
                 x.6h[, sorted.tps+6] + x.6h[, sorted.tps+7] + x.6h[, sorted.tps+8] +
                 x.6h[, sorted.tps+9] + x.6h[, sorted.tps+10] + x.6h[, sorted.tps+11])) # 9-168h
colnames(x.6h) <- c(paste0("before_", marks),
                    paste0("expression_", marks),
                    paste0("after_", marks))

# 12.3. genes expressed at 9 hours
x.9h <- x[x$expression_tp == 3, ]
x.9h <- cbind((x.9h[, sorted.tps] + x.9h[, sorted.tps+1] + x.9h[, sorted.tps+2]), #0-6h
              x.9h[, sorted.tps+3], #9h
              (x.9h[, sorted.tps+4] + x.9h[, sorted.tps+5] + x.9h[, sorted.tps+6] + 
                 x.9h[, sorted.tps+7] + x.9h[, sorted.tps+8] + x.9h[, sorted.tps+9] + 
                 x.9h[, sorted.tps+10] + x.9h[, sorted.tps+11])) # 12-168h
colnames(x.9h) <- c(paste0("before_", marks),
                    paste0("expression_", marks),
                    paste0("after_", marks))

# 12.4. genes expressed at 12 hours
x.12h <- x[x$expression_tp == 4, ]
x.12h <- cbind((x.12h[, sorted.tps] + x.12h[, sorted.tps+1] + x.12h[, sorted.tps+2] + x.12h[, sorted.tps+3]), #0-9h
               x.12h[, sorted.tps+4], #12h
               (x.12h[, sorted.tps+5] + x.12h[, sorted.tps+6] + x.12h[, sorted.tps+7] + 
                  x.12h[, sorted.tps+8] + x.12h[, sorted.tps+9] + x.12h[, sorted.tps+10] + x.12h[, sorted.tps+11])) #18-168h
colnames(x.12h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 12.5. genes expressed at 18 hours
x.18h <- x[x$expression_tp == 5, ]
x.18h <- cbind((x.18h[, sorted.tps] + x.18h[, sorted.tps+1] + x.18h[, sorted.tps+2] + 
                  x.18h[, sorted.tps+3] + x.18h[, sorted.tps+4]), #0-12h
               x.18h[, sorted.tps+5], #18h
               (x.18h[, sorted.tps+6] + x.18h[, sorted.tps+7] + x.18h[, sorted.tps+8] +
                  x.18h[, sorted.tps+9] + x.18h[, sorted.tps+10] + x.18h[, sorted.tps+11])) #24-168h
colnames(x.18h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 12.6. genes expressed at 24 hours
x.24h <- x[x$expression_tp == 6, ]
x.24h <- cbind((x.24h[, sorted.tps] + x.24h[, sorted.tps+1] + x.24h[, sorted.tps+2] + 
                  x.24h[, sorted.tps+3] + x.24h[, sorted.tps+4] + x.24h[, sorted.tps+5]), #0-18h
               x.24h[, sorted.tps+6], #24h
               (x.24h[, sorted.tps+7] + x.24h[, sorted.tps+8] + x.24h[, sorted.tps+9] + 
                  x.24h[, sorted.tps+10] + x.24h[, sorted.tps+11])) #36-168h
colnames(x.24h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 12.7. genes expressed at 36 hours
x.36h <- x[x$expression_tp == 7, ]
x.36h <- cbind((x.36h[, sorted.tps] + x.36h[, sorted.tps+1] + x.36h[, sorted.tps+2] + 
                  x.36h[, sorted.tps+3] + x.36h[, sorted.tps+4] + x.36h[, sorted.tps+5] +
                  x.36h[, sorted.tps+6]), #0-24h
               x.36h[, sorted.tps+7], #36h
               (x.36h[, sorted.tps+8] + x.36h[, sorted.tps+9] + x.36h[, sorted.tps+10] + x.36h[, sorted.tps+11])) #48-168h
colnames(x.36h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 12.8. genes expressed at 48 hours
x.48h <- x[x$expression_tp == 8, ]
x.48h <- cbind((x.48h[, sorted.tps] + x.48h[, sorted.tps+1] + x.48h[, sorted.tps+2] + 
                  x.48h[, sorted.tps+3] + x.48h[, sorted.tps+4] + x.48h[, sorted.tps+5] + 
                  x.48h[, sorted.tps+6] + x.48h[, sorted.tps+7]), #0-36h
               x.48h[, sorted.tps+8], #48h
               (x.48h[, sorted.tps+9] + x.48h[, sorted.tps+10] + x.48h[, sorted.tps+11])) #72-168h
colnames(x.48h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 12.9. genes expressed at 72 hours
x.72h <- x[x$expression_tp == 9, ]
x.72h <- cbind((x.72h[, sorted.tps] + x.72h[, sorted.tps+1] + x.72h[, sorted.tps+2] + 
                  x.72h[, sorted.tps+3] + x.72h[, sorted.tps+4] + x.72h[, sorted.tps+5] + 
                  x.72h[, sorted.tps+6] + x.72h[, sorted.tps+7] + x.72h[, sorted.tps+8]), #0-48h
               x.72h[, sorted.tps+9], #72h
               (x.72h[, sorted.tps+10] + x.72h[, sorted.tps+11])) #120-168h 
colnames(x.72h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 12.10. genes expressed at 120 hours
x.120h <- x[x$expression_tp == 10, ]
x.120h <- cbind((x.120h[, sorted.tps] + x.120h[, sorted.tps+1] + x.120h[, sorted.tps+2] + 
                   x.120h[, sorted.tps+3] + x.120h[, sorted.tps+4] + x.120h[, sorted.tps+5] + 
                   x.120h[, sorted.tps+6] + x.120h[, sorted.tps+7] + x.120h[, sorted.tps+8] + x.120h[, sorted.tps+9]), #0-72h
                x.120h[, sorted.tps+10], #120h
                x.120h[, sorted.tps+11]) #168h 
colnames(x.120h) <- c(paste0("before_", marks),
                      paste0("expression_", marks),
                      paste0("after_", marks))

# 13. rbind all matrices into a unique dataframe
x <- rbind(x.3h, x.6h, x.9h, x.12h, x.18h, x.24h, x.36h, x.48h, x.72h, x.120h)


# 14. if a gene has already acquired mark in a previous time-point
# then it will have the mark also in the next time-point

x$expression_H3K4me1 <- ifelse(x$before_H3K4me1 > 0, 1, x$expression_H3K4me1)
x$after_H3K4me1 <- ifelse(x$expression_H3K4me1 > 0, 1, x$after_H3K4me1)

x$expression_H3K4me2 <- ifelse(x$before_H3K4me2 > 0, 2, x$expression_H3K4me2)
x$after_H3K4me2 <- ifelse(x$expression_H3K4me2 > 0, 2, x$after_H3K4me2)

x$expression_H3K27ac <- ifelse(x$before_H3K27ac > 0, 3, x$expression_H3K27ac)
x$after_H3K27ac <- ifelse(x$expression_H3K27ac > 0, 3, x$after_H3K27ac)

x$expression_H3K9ac <- ifelse(x$before_H3K9ac > 0, 4, x$expression_H3K9ac)
x$after_H3K9ac <- ifelse(x$expression_H3K9ac > 0, 4, x$after_H3K9ac)

x$expression_H3K4me3 <- ifelse(x$before_H3K4me3 > 0, 5, x$expression_H3K4me3)
x$after_H3K4me3 <- ifelse(x$expression_H3K4me3 > 0, 5, x$after_H3K4me3)

x$expression_H3K36me3 <- ifelse(x$before_H3K36me3 > 0, 6, x$expression_H3K36me3)
x$after_H3K36me3 <- ifelse(x$expression_H3K36me3 > 0, 6, x$after_H3K36me3)

x$expression_H4K20me1 <- ifelse(x$before_H4K20me1 > 0, 7, x$expression_H4K20me1)
x$after_H4K20me1 <- ifelse(x$expression_H4K20me1 > 0, 7, x$after_H4K20me1)

x$expression_H3K9me3 <- ifelse(x$before_H3K9me3 > 0, 8, x$expression_H3K9me3)
x$after_H3K9me3 <- ifelse(x$expression_H3K9me3 > 0, 8, x$after_H3K9me3)

x$expression_H3K27me3 <- ifelse(x$before_H3K27me3 > 0, 9, x$expression_H3K27me3)
x$after_H3K27me3 <- ifelse(x$expression_H3K27me3 > 0, 9, x$after_H3K27me3)


# 15. genes unmarked by H3K27me3
x1.ser.proxEnhancers <- x[rownames(x) %in% rownames(x1.ser.promoters), ]
# ser.obj1 <- seriate(as.matrix(x1.ser.proxEnhancers[, c(1, 3, 10, 12, 19, 21)]), method = "PCA")
ser.obj1 <- seriate(as.matrix(x1.ser.proxEnhancers[, c(1:5, 10:14, 19:23)]), method = "PCA")
x1.ser.proxEnhancers <- x1.ser.proxEnhancers[get_order(ser.obj1, 1), ]


# 16. genes marked by H3K27me3
x2.ser.proxEnhancers <- x[rownames(x) %in% rownames(x2.ser.promoters), ]
# ser.obj2 <- seriate(as.matrix(x2.ser.proxEnhancers[, c(1, 3, 10, 12, 19, 21)]), method = "PCA")
ser.obj2 <- seriate(as.matrix(x2.ser.proxEnhancers[, c(1:5, 10:14, 19:23)]), method = "PCA")
x2.ser.proxEnhancers <- x2.ser.proxEnhancers[get_order(ser.obj2, 1), ]

x.proxEnhancers <- x
x2.proxEnhancers <- x2.ser.proxEnhancers

# Distal enhancers ----
# 1. set working directory
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
genes.pairs <- read.delim("/users/project/encode_005982_no_backup/flagship/Borsari_et_al/analysis/enhancers/ENCODE-rE2G.blood/interacting_pairs_filtered_100bp.merged_closest.bed",
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


# 15. add to x the info on gene-expression time-point
stopifnot(identical(rownames(x), rownames(expression.tp)))
x$expression_tp <- expression.tp$expression_tp


# 16. columns corresponding to the seven active marks
sorted.tps <- c(1, # H3K4me1
                13, # H3K4me2
                25, # H3K27ac
                37, # H3K9ac
                49, # H3K4me3
                61, # H3K36me3
                73, # H4K20me1
                85, # H3K9me3
                97) # H3K27me3
#marks <- marks[-9]


#----------------------------------------------------
# 17. subset genes based on time-point of expression
#----------------------------------------------------

# 17.1. genes expressed at 3 hours
x.3h <- x[x$expression_tp == 1, ]
x.3h <- cbind(x.3h[, sorted.tps], # 0h 
              x.3h[, sorted.tps+1], # 3h
              (x.3h[, sorted.tps+2] + x.3h[, sorted.tps+3] + x.3h[, sorted.tps+4] + 
                 x.3h[, sorted.tps+5] + x.3h[, sorted.tps+6] + x.3h[, sorted.tps+7] + 
                 x.3h[, sorted.tps+8] + x.3h[, sorted.tps+9] + x.3h[, sorted.tps+10] + x.3h[, sorted.tps+11])) # 6-168h
colnames(x.3h) <- c(paste0("before_", marks),
                    paste0("expression_", marks),
                    paste0("after_", marks))

# 17.2. genes expressed at 6 hours
x.6h <- x[x$expression_tp == 2, ]
x.6h <- cbind((x.6h[, sorted.tps] + x.6h[, sorted.tps+1]), # 0-3h 
              x.6h[, sorted.tps+2], # 6h
              (x.6h[, sorted.tps+3] + x.6h[, sorted.tps+4] + x.6h[, sorted.tps+5] + 
                 x.6h[, sorted.tps+6] + x.6h[, sorted.tps+7] + x.6h[, sorted.tps+8] +
                 x.6h[, sorted.tps+9] + x.6h[, sorted.tps+10] + x.6h[, sorted.tps+11])) # 9-168h
colnames(x.6h) <- c(paste0("before_", marks),
                    paste0("expression_", marks),
                    paste0("after_", marks))

# 17.3. genes expressed at 9 hours
x.9h <- x[x$expression_tp == 3, ]
x.9h <- cbind((x.9h[, sorted.tps] + x.9h[, sorted.tps+1] + x.9h[, sorted.tps+2]), #0-6h
              x.9h[, sorted.tps+3], #9h
              (x.9h[, sorted.tps+4] + x.9h[, sorted.tps+5] + x.9h[, sorted.tps+6] + 
                 x.9h[, sorted.tps+7] + x.9h[, sorted.tps+8] + x.9h[, sorted.tps+9] + 
                 x.9h[, sorted.tps+10] + x.9h[, sorted.tps+11])) # 12-168h
colnames(x.9h) <- c(paste0("before_", marks),
                    paste0("expression_", marks),
                    paste0("after_", marks))

# 17.4. genes expressed at 12 hours
x.12h <- x[x$expression_tp == 4, ]
x.12h <- cbind((x.12h[, sorted.tps] + x.12h[, sorted.tps+1] + x.12h[, sorted.tps+2] + x.12h[, sorted.tps+3]), #0-9h
               x.12h[, sorted.tps+4], #12h
               (x.12h[, sorted.tps+5] + x.12h[, sorted.tps+6] + x.12h[, sorted.tps+7] + 
                  x.12h[, sorted.tps+8] + x.12h[, sorted.tps+9] + x.12h[, sorted.tps+10] + x.12h[, sorted.tps+11])) #18-168h
colnames(x.12h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 17.5. genes expressed at 18 hours
x.18h <- x[x$expression_tp == 5, ]
x.18h <- cbind((x.18h[, sorted.tps] + x.18h[, sorted.tps+1] + x.18h[, sorted.tps+2] + 
                  x.18h[, sorted.tps+3] + x.18h[, sorted.tps+4]), #0-12h
               x.18h[, sorted.tps+5], #18h
               (x.18h[, sorted.tps+6] + x.18h[, sorted.tps+7] + x.18h[, sorted.tps+8] +
                  x.18h[, sorted.tps+9] + x.18h[, sorted.tps+10] + x.18h[, sorted.tps+11])) #24-168h
colnames(x.18h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 17.6. genes expressed at 24 hours
x.24h <- x[x$expression_tp == 6, ]
x.24h <- cbind((x.24h[, sorted.tps] + x.24h[, sorted.tps+1] + x.24h[, sorted.tps+2] + 
                  x.24h[, sorted.tps+3] + x.24h[, sorted.tps+4] + x.24h[, sorted.tps+5]), #0-18h
               x.24h[, sorted.tps+6], #24h
               (x.24h[, sorted.tps+7] + x.24h[, sorted.tps+8] + x.24h[, sorted.tps+9] + 
                  x.24h[, sorted.tps+10] + x.24h[, sorted.tps+11])) #36-168h
colnames(x.24h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 17.7. genes expressed at 36 hours
x.36h <- x[x$expression_tp == 7, ]
x.36h <- cbind((x.36h[, sorted.tps] + x.36h[, sorted.tps+1] + x.36h[, sorted.tps+2] + 
                  x.36h[, sorted.tps+3] + x.36h[, sorted.tps+4] + x.36h[, sorted.tps+5] +
                  x.36h[, sorted.tps+6]), #0-24h
               x.36h[, sorted.tps+7], #36h
               (x.36h[, sorted.tps+8] + x.36h[, sorted.tps+9] + x.36h[, sorted.tps+10] + x.36h[, sorted.tps+11])) #48-168h
colnames(x.36h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 17.8. genes expressed at 48 hours
x.48h <- x[x$expression_tp == 8, ]
x.48h <- cbind((x.48h[, sorted.tps] + x.48h[, sorted.tps+1] + x.48h[, sorted.tps+2] + 
                  x.48h[, sorted.tps+3] + x.48h[, sorted.tps+4] + x.48h[, sorted.tps+5] + 
                  x.48h[, sorted.tps+6] + x.48h[, sorted.tps+7]), #0-36h
               x.48h[, sorted.tps+8], #48h
               (x.48h[, sorted.tps+9] + x.48h[, sorted.tps+10] + x.48h[, sorted.tps+11])) #72-168h
colnames(x.48h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 17.9. genes expressed at 72 hours
x.72h <- x[x$expression_tp == 9, ]
x.72h <- cbind((x.72h[, sorted.tps] + x.72h[, sorted.tps+1] + x.72h[, sorted.tps+2] + 
                  x.72h[, sorted.tps+3] + x.72h[, sorted.tps+4] + x.72h[, sorted.tps+5] + 
                  x.72h[, sorted.tps+6] + x.72h[, sorted.tps+7] + x.72h[, sorted.tps+8]), #0-48h
               x.72h[, sorted.tps+9], #72h
               (x.72h[, sorted.tps+10] + x.72h[, sorted.tps+11])) #120-168h 
colnames(x.72h) <- c(paste0("before_", marks),
                     paste0("expression_", marks),
                     paste0("after_", marks))


# 17.10. genes expressed at 120 hours
x.120h <- x[x$expression_tp == 10, ]
x.120h <- cbind((x.120h[, sorted.tps] + x.120h[, sorted.tps+1] + x.120h[, sorted.tps+2] + 
                   x.120h[, sorted.tps+3] + x.120h[, sorted.tps+4] + x.120h[, sorted.tps+5] + 
                   x.120h[, sorted.tps+6] + x.120h[, sorted.tps+7] + x.120h[, sorted.tps+8] + x.120h[, sorted.tps+9]), #0-72h
                x.120h[, sorted.tps+10], #120h
                x.120h[, sorted.tps+11]) #168h 
colnames(x.120h) <- c(paste0("before_", marks),
                      paste0("expression_", marks),
                      paste0("after_", marks))

# 18. rbind all matrices into a unique dataframe
x <- rbind(x.3h, x.6h, x.9h, x.12h, x.18h, x.24h, x.36h, x.48h, x.72h, x.120h)


# 19. if a gene has already acquired a mark in a previous time-point
# then it will have the mark also in the next time-points
x$expression_H3K4me1 <- ifelse(x$before_H3K4me1 > 0, 1, x$expression_H3K4me1)
x$after_H3K4me1 <- ifelse(x$expression_H3K4me1 > 0, 1, x$after_H3K4me1)

x$expression_H3K4me2 <- ifelse(x$before_H3K4me2 > 0, 2, x$expression_H3K4me2)
x$after_H3K4me2 <- ifelse(x$expression_H3K4me2 > 0, 2, x$after_H3K4me2)

x$expression_H3K27ac <- ifelse(x$before_H3K27ac > 0, 3, x$expression_H3K27ac)
x$after_H3K27ac <- ifelse(x$expression_H3K27ac > 0, 3, x$after_H3K27ac)

x$expression_H3K9ac <- ifelse(x$before_H3K9ac > 0, 4, x$expression_H3K9ac)
x$after_H3K9ac <- ifelse(x$expression_H3K9ac > 0, 4, x$after_H3K9ac)

x$expression_H3K4me3 <- ifelse(x$before_H3K4me3 > 0, 5, x$expression_H3K4me3)
x$after_H3K4me3 <- ifelse(x$expression_H3K4me3 > 0, 5, x$after_H3K4me3)

x$expression_H3K36me3 <- ifelse(x$before_H3K36me3 > 0, 6, x$expression_H3K36me3)
x$after_H3K36me3 <- ifelse(x$expression_H3K36me3 > 0, 6, x$after_H3K36me3)

x$expression_H4K20me1 <- ifelse(x$before_H4K20me1 > 0, 7, x$expression_H4K20me1)
x$after_H4K20me1 <- ifelse(x$expression_H4K20me1 > 0, 7, x$after_H4K20me1)

x$expression_H3K9me3 <- ifelse(x$before_H3K9me3 > 0, 8, x$expression_H3K9me3)
x$after_H3K9me3 <- ifelse(x$expression_H3K9me3 > 0, 8, x$after_H3K9me3)

x$expression_H3K27me3 <- ifelse(x$before_H3K27me3 > 0, 9, x$expression_H3K27me3)
x$after_H3K27me3 <- ifelse(x$expression_H3K27me3 > 0, 9, x$after_H3K27me3)


# 20. genes unmarked by H3K27me3
x1.ser.distEnhancers <- x[rownames(x) %in% rownames(x1.ser.promoters), ]
# x1.ser.distEnhancers <- x1.ser.distEnhancers[complete.cases(x1.ser.distEnhancers), ]
x1.ser.distEnhancers[is.na(x1.ser.distEnhancers)] <- -1
# ser.obj1 <- seriate(as.matrix(x1.ser.distEnhancers[, c(1, 3, 10, 12, 19, 21)]), method = "PCA")
ser.obj1 <- seriate(as.matrix(x1.ser.distEnhancers[, c(1:5, 10:14, 19:23)]), method = "PCA")
x1.ser.distEnhancers <- x1.ser.distEnhancers[get_order(ser.obj1, 1), ]


# 21. genes marked by H3K27me3
x2.ser.distEnhancers <- x[rownames(x) %in% rownames(x2.ser.promoters), ]
# x2.ser.distEnhancers <- x2.ser.distEnhancers[complete.cases(x2.ser.distEnhancers), ]
x2.ser.distEnhancers[is.na(x2.ser.distEnhancers)] <- -1
# ser.obj2 <- seriate(as.matrix(x2.ser.distEnhancers[, c(1, 3, 10, 12, 19, 21)]), method = "PCA")
ser.obj2 <- seriate(as.matrix(x2.ser.distEnhancers[, c(1:5, 10:14, 19:23)]), method = "PCA")
x2.ser.distEnhancers <- x2.ser.distEnhancers[get_order(ser.obj2, 1), ]

x.distEnhancers <- x
x2.distEnhancers <- x2.ser.distEnhancers

# Integrate data from promoters, proximal and distal enhancers ----
# source("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/paper.figures/figure.5b.pEnh.R")
rm(list=setdiff(ls(), c("x1.ser.proxEnhancers", "x2.ser.proxEnhancers",
                        "x1.ser.promoters", "x2.ser.promoters",
                        "x1.ser.proxEnhancers", "x2.ser.proxEnhancers",
                        "x1.ser.distEnhancers", "x2.ser.distEnhancers")))


# Sort all matrices based on the same order - new addition 09/01/2026
## Generate matrix binding all 3 original matrices
tmp <- do.call("cbind", list(x1.ser.promoters, 
                             x1.ser.proxEnhancers[rownames(x1.ser.promoters), ], 
                             x1.ser.distEnhancers[rownames(x1.ser.promoters), ]))

## Convert NAs to 0 for genes without known distal enhancer
tmp[is.na(tmp)] <- 0
ser.obj1 <- seriate(as.matrix(tmp[, c(1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 19, 21, 22, 23, 24,
                                      28, 29, 30, 31, 32, 37, 38, 39, 40, 41, 46, 47, 48, 49, 50,
                                      55, 56, 57, 58, 59, 64, 65, 66, 67, 68, 73, 74, 75, 76, 77)]), method = "PCA")
                                                                                                             
x1.ser.promoters <- x1.ser.promoters[get_order(ser.obj1, 1), ]
x1.ser.proxEnhancers <- x1.ser.proxEnhancers[rownames(x1.ser.promoters[get_order(ser.obj1, 1), ]), ]
x1.ser.distEnhancers <- x1.ser.distEnhancers[rownames(x1.ser.promoters[get_order(ser.obj1, 1), ]), ]

## Generate matrix binding all 3 original matrices
tmp <- do.call("cbind", list(x2.ser.promoters, 
                             x2.ser.proxEnhancers[rownames(x2.ser.promoters), ], 
                             x2.ser.distEnhancers[rownames(x2.ser.promoters), ]))

## Convert NAs to 0 for genes without known distal enhancer
tmp[is.na(tmp)] <- 0
ser.obj2 <- seriate(as.matrix(tmp[, c(1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 19, 21, 22, 23, 24,
                                      28, 29, 30, 31, 32, 37, 38, 39, 40, 41, 46, 47, 48, 49, 50,
                                      55, 56, 57, 58, 59, 64, 65, 66, 67, 68, 73, 74, 75, 76, 77)]), method = "PCA")

x2.ser.promoters <- x2.ser.promoters[get_order(ser.obj2, 1), ]
x2.ser.proxEnhancers <- x2.ser.proxEnhancers[rownames(x2.ser.promoters[get_order(ser.obj2, 1), ]), ]
x2.ser.distEnhancers <- x2.ser.distEnhancers[rownames(x2.ser.promoters[get_order(ser.obj2, 1), ]), ]


pdf("fig_4f.pdf",
    height = 2.5,
    width = 4)

#--------------------
# promoters 
#--------------------

# before expression
 h <- Heatmap(rbind(x1.ser.promoters[1:145, c(1:9)], x2.ser.promoters[1:112, c(1:9)]), 
        column_title = "before\nexpression",
        column_title_gp = gpar(fontsize = 10),
        name = "promoter 1", 
        cluster_rows = F,
        cluster_columns = F,
        col= c("white", "#E5AB00", "#A67C00", "#630039", "#AF4D85", "#D199B9",
               "#7FBC41", "#4C7027", "#A7ADD4", "#1D2976"),
        show_heatmap_legend = F,
        show_row_names = FALSE, 
        width = unit(25, "mm"),
        show_column_names = F,
        gap = unit(1, "mm"),
        split = c(rep("a", 145), 
                  rep("b", 112)),
        border = T,
        use_raster = TRUE) +
  
  # with expression
  Heatmap(rbind(x1.ser.promoters[1:145, c(10:18)], x2.ser.promoters[1:112, c(10:18)]),
          column_title = "with\nexpression",
          column_title_gp = gpar(fontsize = 10),
          name = "promoter 2", 
          cluster_rows = F,
          cluster_columns = F,
          col= c("white", "#E5AB00", "#A67C00", "#630039", "#AF4D85", "#D199B9",
                 "#7FBC41", "#4C7027", "#A7ADD4", "#1D2976"),
          show_heatmap_legend = F,
          show_row_names = FALSE, 
          width = unit(25, "mm"),
          show_column_names = F,
          gap = unit(1, "mm"),
          split = c(rep("a", 145), 
                    rep("b", 112)),
          border = T,
          use_raster = TRUE) +
  
  # after expression
  Heatmap(rbind(x1.ser.promoters[1:145, c(19:27)], x2.ser.promoters[1:112, c(19:27)]),
          column_title = "after\nexpression",
          column_title_gp = gpar(fontsize = 10),
          name = "promoter 3", 
          cluster_rows = F,
          cluster_columns = F,
          col= c("white", "#E5AB00", "#A67C00", "#630039", "#AF4D85", "#D199B9",
                 "#7FBC41", "#4C7027", "#A7ADD4", "#1D2976"),
          show_heatmap_legend = F,
          show_row_names = FALSE, 
          width = unit(25, "mm"),
          show_column_names = F,
          gap = unit(1, "mm"),
          split = c(rep("a", 145), 
                    rep("b", 112)),
          border = T,
          use_raster = TRUE)

ComplexHeatmap::draw(h, row_title = "promoters")

#----------
# proximal enhancers
#-----------

# before expression
h <- Heatmap(rbind(x1.ser.proxEnhancers[145:1, c(1:9)], x2.ser.proxEnhancers[1:112, c(1:9)]), 
        column_title = "before\nexpression",
        column_title_gp = gpar(fontsize = 10),
        name = "proximal region 1", 
        cluster_rows = F,
        cluster_columns = F,
        col= c("white", "#E5AB00", "#A67C00", "#630039", "#AF4D85", "#D199B9",
               "#7FBC41", "#4C7027", "#A7ADD4", "#1D2976"),
        show_heatmap_legend = F,
        show_row_names = FALSE, 
        width = unit(25, "mm"),
        show_column_names = F,
        gap = unit(1, "mm"),
        split = c(rep("a", 145), 
                  rep("b", 112)),
        border = T,
        use_raster = TRUE) +
  
  # with expression
  Heatmap(rbind(x1.ser.proxEnhancers[145:1, c(10:18)], x2.ser.proxEnhancers[1:112, c(10:18)]),
          column_title = "with\nexpression",
          column_title_gp = gpar(fontsize = 10),
          name = "proximal region 2", 
          cluster_rows = F,
          cluster_columns = F,
          col= c("white", "#E5AB00", "#A67C00", "#630039", "#AF4D85", "#D199B9",
                 "#7FBC41", "#4C7027", "#A7ADD4", "#1D2976"),
          show_heatmap_legend = F,
          show_row_names = FALSE, 
          width = unit(25, "mm"),
          show_column_names = F,
          gap = unit(1, "mm"),
          split = c(rep("a", 145), 
                    rep("b", 112)),
          border = T,
          use_raster = TRUE) +
  
  # after expression
  Heatmap(rbind(x1.ser.proxEnhancers[145:1, c(19:27)], x2.ser.proxEnhancers[1:112, c(19:27)]),
          column_title = "after\nexpression",
          column_title_gp = gpar(fontsize = 10),
          name = "proximal region 3", 
          cluster_rows = F,
          cluster_columns = F,
          col= c("white", "#E5AB00", "#A67C00", "#630039", "#AF4D85", "#D199B9",
                 "#7FBC41", "#4C7027", "#A7ADD4", "#1D2976"),
          show_heatmap_legend = F,
          show_row_names = FALSE, 
          width = unit(25, "mm"),
          show_column_names = F,
          gap = unit(1, "mm"),
          split = c(rep("a", 145), 
                    rep("b", 112)),
          border = T,
          use_raster = TRUE)

ComplexHeatmap::draw(h, row_title = "proximal enhancers")

#-----------
# distal enhancers
#-----------

# before expression
h <- Heatmap(rbind(x1.ser.distEnhancers[145:1, c(1:9)], x2.ser.distEnhancers[1:112, c(1:9)]), 
        column_title = "before\nexpression",
        column_title_gp = gpar(fontsize = 10),
        name = "distal enhancer 1", 
        cluster_rows = F,
        cluster_columns = F,
        na_col = "white", 
        col= c("-1" = "grey", "0" = "white", "1" = "#E5AB00", 
               "2" = "#A67C00", "3" = "#630039", "4" = "#AF4D85", "5" = "#D199B9",
               "6" = "#7FBC41", "7" = "#4C7027", "8" = "#A7ADD4", "9" = "#1D2976"),
        show_heatmap_legend = F,
        show_row_names = FALSE, 
        width = unit(25, "mm"),
        show_column_names = F,
        gap = unit(1, "mm"),
        # split = c(rep("a", 133),
        #           rep("b", 108)),
        split = c(rep("a", 145),
                  rep("b", 112)),
        border = TRUE,
        use_raster = TRUE) +
  
  # with expression
  Heatmap(rbind(x1.ser.distEnhancers[145:1, c(10:18)], x2.ser.distEnhancers[1:112, c(10:18)]),
          column_title = "with\nexpression",
          column_title_gp = gpar(fontsize = 10),
          name = "distal enhancer 2", 
          cluster_rows = F,
          cluster_columns = F,
          na_col = "white", 
          col= c("-1" = "grey", "0" = "white", "1" = "#E5AB00", 
                 "2" = "#A67C00", "3" = "#630039", "4" = "#AF4D85", "5" = "#D199B9",
                 "6" = "#7FBC41", "7" = "#4C7027", "8" = "#A7ADD4", "9" = "#1D2976"),
          show_heatmap_legend = F,
          show_row_names = FALSE, 
          width = unit(25, "mm"),
          show_column_names = F,
          gap = unit(1, "mm"),
          # split = c(rep("a", 133),
          #           rep("b", 108)),
          split = c(rep("a", 145),
                    rep("b", 112)),
          border = TRUE,
          use_raster = TRUE) +
  
  # after expression
  Heatmap(rbind(x1.ser.distEnhancers[145:1, c(19:27)], x2.ser.distEnhancers[1:112, c(19:27)]), 
          column_title = "after\nexpression",
          column_title_gp = gpar(fontsize = 10),
          name = "distal enhancer 3", 
          cluster_rows = F,
          cluster_columns = F,
          na_col = "white", 
          col= c("-1" = "grey", "0" = "white", "1" = "#E5AB00", 
                 "2" = "#A67C00", "3" = "#630039", "4" = "#AF4D85", "5" = "#D199B9",
                 "6" = "#7FBC41", "7" = "#4C7027", "8" = "#A7ADD4", "9" = "#1D2976"),
          show_heatmap_legend = F,
          show_row_names = FALSE, 
          width = unit(25, "mm"),
          show_column_names = F,
          gap = unit(1, "mm"),
          # split = c(rep("a", 133),
          #           rep("b", 108)),
          split = c(rep("a", 145),
                    rep("b", 112)),
          border = TRUE,
          use_raster = TRUE)

ComplexHeatmap::draw(h, row_title = "distal enhancers")

dev.off()
