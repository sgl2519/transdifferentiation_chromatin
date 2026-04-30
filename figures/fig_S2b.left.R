library(pheatmap)
library(RColorBrewer)

# Set function for finding consecutive elements
find_consecutive_with_values <- function(vec) {
  # Initialize result vectors for counts and values
  counts <- numeric(0)
  values <- numeric(0)
  
  # Start with the first element
  current_value <- vec[1]
  count <- 1
  
  # Loop through the vector starting from the second element
  for (i in 2:length(vec)) {
    if (vec[i] == current_value) {
      # If the same as the previous one, increment the count
      count <- count + 1
    } else {
      # If different, append the count and the value to the results
      counts <- c(counts, count)
      values <- c(values, current_value)
      count <- 1
      current_value <- vec[i]
    }
  }
  
  # Append the last group count and value
  counts <- c(counts, count)
  values <- c(values, current_value)
  
  return(list(counts = counts, values = values))
}


# Load cCREs that gain marking
setwd("../analysis/cCREs/precedence/")

marks <- c("H3K27ac", "H3K9ac", "H3K4me3", 
           "H3K4me1", "H3K4me2", "H3K36me3", 
           "H4K20me1", "H3K9me3", "H3K27me3", "ATAC-seq")

# Load the matrix of mark appearance
m <- read.table('mark.gain_tp.tsv',
                header = T)

# Remove the cases where no active mark appears during transdifferentiation
m <- m[rowSums(m[, c(3:9)]) > 0, ]

# Remove the cases where no active mark appears during transdifferentiation (without considering H3K4me1/2)
#m <- m[rowSums(m[, c(3:5, 8:9)]) > 0, ]
m <- m[rowSums(m[, c(3:5)]) > 0, ]

# Select only those cases where all active marks appear (i.e, are not already present) during the process
#m <- m[rowSums(m[, c(3:5, 8:9)] == 1)==0, , drop = FALSE]
m <- m[rowSums(m[, c(3:5)] == 1)==0, , drop = FALSE]
m <- m[rowSums(m[, c(3:5)] == 2)==0, , drop = FALSE]

# Select cases that when gaining marking, they keep it for a determined (cons_tps) number of time points
cons_tps <- 3

consistent <- list()
for ( i in marks[1:5] ) {
  # Read peak presence and absence matrix
  tmp <- read.table(paste0(i,
                           '/',
                           i,
                           '.peaks.dynamics.binary.tsv'))
  
  # Subset the cCREs from the previous matrix
  tmp <- tmp[tmp$V1 %in% m$cCRE_id, ]
  
  
  # Calculate the number of consecutive time points with marking
  tmp$consecutive <- apply(tmp[, 2:13], 1, function(a) {
    x <- find_consecutive_with_values(as.integer(a))
    if( any(x$counts[x$values == 1] >= cons_tps) ) {
      return('consistent')
    } else {
      return('inconsistent')
    }
    
  })
  
  # Save the cCRE ids corresponding to the consistent marking by the corresponding mark
  consistent[[i]] <- tmp[tmp$consecutive == 'consistent', ]$V1
  
}

# Subset the cases that gain a mark consistently for the set of activating marks
m <- m[m$cCRE_id %in% Reduce(union, list(consistent$H3K27ac, 
                                         consistent$H3K9ac, 
                                         consistent$H3K4me3)), ]

# Add ATAC-seq data
c <- read.table('ATACseq/cCREs.ATACseq_overlap.tsv')

## Convert it to short format
c <- c[, c(6:8, 13)]
c <- reshape(c, idvar = c('V6', 'V7'), timevar = "V8", direction = "wide")
c <- c[, -3]
c <- c[, c('V6', 'V7', 'V13.H000', 'V13.H012', 'V13.H024', 'V13.H096')]
c[is.na(c)] <- 0
rownames(c) <- c$V7
c <- c[, 3:6]
c[c > 0] <- 1 

## Obtain time point of first appearance
### Determine the equivalente between timepoints in ChIP-seq and ATAC-seq
tps <- c(1, 5, 7, 10.5)
c$appear <- apply(c, 1, function(a) tps[which(a!=0)[1]])
c[is.na(c)] <- 0
c$cCRE_id <- rownames(c)

n <- m

n[is.na(n)] <- 0
n$Row.names <- NULL

## Obtain time point of first appearance
tps <- 1:12
n$appear <- apply(n[, 2:12], 1, function(a) tps[which(a!=0)[1]])
n[is.na(n)] <- 0

## Remove cases that do not gain ATAC-seq peak
c <- c[c$appear != 0, ]

# Merge with ATAC-seq
m <- merge(m, c[, 5:6], by = 'cCRE_id', all = F)
colnames(m) <- c(colnames(m)[1:11], 'ATAC-seq')

# Merge with cEBPa
m <- merge(m, n[, c(1, 14)], by = 'cCRE_id', all = F)
colnames(m) <- c(colnames(m)[1:12], 'cEBPa')

# Divide the dataset according to categories
m$supracategory <- 'CA'
m[m$category == 'PLS', ]$supracategory <- 'PLS'
m[m$category == 'pELS', ]$supracategory <- 'pELS'
m[m$category == 'dELS', ]$supracategory <- 'dELS'

# Substitute 0 with 13
m[m == 0] <- 13

# Reorder data frame columns
m <- m[, c(1:9, 12, 13, 10, 11, 14)]

pdf("fig_S2b.left.pdf",
    width = 5,
    height = 5)
for ( i in unique(m$supracategory) ) {
  # Generate matrix of precedence
  d <- matrix(, nrow = 9, ncol = 9)
  colnames(d) <- c(marks[c(1:7, 10)], 'cEBPa')
  rownames(d) <- c(marks[c(1:7, 10)], 'cEBPa')
  
  # Fill matrix of precedence
  for ( j in c(3:11) ) {
    for ( k in c(3:11) ) {
      v <- m[m$supracategory == i, j] - m[m$supracategory == i, k]
      v <- sign(v)
      
      # Count number of cases of precedence
      p <- sum(v == -1)
      
      # Count total number of cases
      t <- length(v)
      
      # Add to the matrix
      d[j - 2, k - 2] <- p/t*100
      
    }
  }
  
  if (i == 'PLS') {
    order <- c("H3K4me1", "H3K4me2",
               "H3K27ac", "ATAC-seq", 
               "H3K9ac", "H3K4me3",
               "cEBPa",
               "H4K20me1", "H3K36me3")
    
    
  }
  
  if (i == 'pELS') {
    order <- c("H3K4me1", "H3K4me2",
               "H3K27ac", "ATAC-seq", 
               "H3K9ac",
               "cEBPa", "H3K4me3",
               "H4K20me1", "H3K36me3")
    
    
  }
  
  if (i == 'dELS') {
    order <- c("H3K4me1", "H3K4me2",
               "H3K27ac", "ATAC-seq",  
               "H3K9ac",
               "cEBPa", "H3K4me3",
               "H4K20me1", "H3K36me3")
    
    
  }
  
  if (i == 'CA') {
    order <- c("H3K4me1", "H3K4me2",
               "H3K27ac", "ATAC-seq",
               "cEBPa", 
               "H3K9ac", "H3K4me3",
               "H4K20me1", "H3K36me3")
    
    
  }
  
  d <- d[order, order]
  
  pheatmap::pheatmap(d, 
                     main = paste0('Precedence for ', 
                                   nrow(m[m$supracategory == i, ]),
                                   ' ',
                                   i, 's'),
                     cluster_rows = F,
                     cluster_cols = F,
                     display_numbers = T,
                     border_color = NA,
                     fontsize = 15,
                     legend = F, 
                     number_format = "%.1f",
                     color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(19),
                     breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90))
}

dev.off()
