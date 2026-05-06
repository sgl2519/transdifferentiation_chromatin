
#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c("--rep1"), default = NULL, help = "Replicate n. 1. WARNING: colnames of rep. 1 and 2 must be in the same order" ),
  
  make_option( c("--rep2"), default = NULL, help = "Replicate n. 2" ),
  
  make_option( c("-o", "--output"), default = "out.maSigPro.tsv",
               help="Output file name. 'stdout' for standard output [default=%default].")
  
)


parser <- OptionParser(
  usage = "%prog [options] files", 
  option_list=option_list,
  description = "\nWrapper for maSigPro (significant genes detected at 5% FDR)."
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options




#************
# LIBRARIES *
#************

library(maSigPro)
library(dplyr)
library(tidyr)



#**************
# READ OPTIONS
#**************

# debugging options
# a <- read.table("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/H3K4me3/QN.merged/H3K4me3.R1.matrix.after.QN.merged.tsv",
#                  h=T, sep="\t", row.names = 1)
# b <- read.table("/no_backup/rg/bborsari/projects/ERC/human/2018-01-19.chip-nf/Borsari_et_al/analysis/all.marks/H3K4me3/QN.merged/H3K4me3.R2.matrix.after.QN.merged.tsv",
#                 h=T, sep="\t", row.names = 1)


if ( ! is.null(opt$rep1) && ! is.null(opt$rep2) ) {
  
  a <- read.table( file = opt$rep1, header = T, quote = NULL, row.names = 1 )
  b <- read.table( file = opt$rep2, header = T, quote = NULL, row.names = 1 )
  
} else {
  
  print('Missing input matrices!')
  quit(save = 'no')
  
}

a$gene_id <- rownames(a)
b$gene_id <- rownames(b)


# output
output = ifelse(opt$output == "stdout", "", opt$output)


#********
# BEGIN *
#********

# merge rep. 1 and rep. 2
a.b.merged <- merge(a, b, by = "gene_id")
rownames(a.b.merged) <- a.b.merged$gene_id
a.b.merged$gene_id <- NULL

# add pseudo-count of 10^-16 for genes
# that have 0 signal
a.b.merged$sum <- apply(a.b.merged, 1, sum)
a.b.merged[a.b.merged$sum == 0, ] <- a.b.merged[a.b.merged$sum == 0, ] + 10^-16
a.b.merged$sum <- NULL

# generate design
Time <- rep(c(0,3,6,9,12,18,24,36,48,72,120,168),2)
Replicates <- rep(c(1,2), each = 12)
Group <- rep(1,24)
ss.edesign <- cbind(Time,Replicates,Group)
rownames(ss.edesign) <- colnames(a.b.merged)
design <- make.design.matrix(ss.edesign, degree = 2)

# fit 
fit <- p.vector(a.b.merged, design, Q = 0.05, MT.adjust = "BH", min.obs = 20, counts = FALSE)
significant.genes <- data.frame(gene_id = rownames(a.b.merged), 
                                p_adjusted = fit$p.adjusted)
rownames(significant.genes) <- significant.genes$gene_id
significant.genes <- significant.genes[significant.genes$p_adjusted < 0.05, "p_adjusted", drop=F]




#*********
# OUTPUT *
#*********

write.table(significant.genes, file=output, quote=FALSE, sep="\t", row.names = T, col.names = T)



#*******
# EXIT *
#*******

quit(save="no")
