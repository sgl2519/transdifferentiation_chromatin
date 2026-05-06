#*****************
# OPTION PARSING *
#*****************

suppressPackageStartupMessages(library("optparse"))

option_list <- list (
  
  make_option( c("-a", "--matrix_A"), default = NULL, help = "Matrix n. 1" ),
  
  make_option( c("-b", "--matrix_B"), default = NULL, help = "Matrix n. 2" ),
  
  make_option( c("--header"), default = TRUE,
               help = "Whether the input matrix has a header [default = %default]." ),
  
  make_option( c("--log10", "-l"), default = FALSE,
               help = "Apply log10 to the matrix as pre-processing step [default = %default]." ),
  
  make_option( c("--pseudocount", "-p"), type = 'numeric', default = 0.001,
               help = "The pseudocount to add when applying the log [default = %default]." ),
  
  make_option( c("--output", "-o"), default = "mean.matrix.tsv",
               help="Output file name. 'stdout' for standard output [default=%default]."),
  
  make_option( c("--round", "-r"), type="numeric", default = NULL,
               help="The number of decimal digits in the final matrix [default=%default]."),
  
  make_option( c("--remove_NAs"), default = TRUE,
               help = "Whether to ignore NAs; if FALSE, the mean between a number and NA will be NA [default=%default].")
  
)


parser <- OptionParser(
  usage = "%prog [options] files", 
  option_list=option_list,
  description = "\nComputes the mean over the 2 input matrices"
)

arguments <- parse_args(parser, positional_arguments = TRUE)
opt <- arguments$options





#***************
# READ OPTIONS *
#***************


 if ( ! is.null(opt$matrix_A) && ! is.null(opt$matrix_B) ) {
  
  a <- read.table( file = opt$matrix_A, header = opt$header, quote = NULL )
  b <- read.table( file = opt$matrix_B, header = opt$header, quote = NULL )
  
} else {
  
  print('Missing input matrices!')
  quit(save = 'no')
  
}


if ( opt$header == FALSE ) {
  
  rownames(a) <- a$V1
  rownames(b) <- b$V1
  a$V1 <- NULL
  b$V1 <- NULL
  
}


a <- as.matrix(a)
b <- as.matrix(b)

if ( !is.numeric(a) && !(is.numeric(b)) ) {
  
  print ("input matrices must contain numeric values")
  quit (save = 'no')
  
}

if ( opt$log10 ) {
  
  a <- log10 ( a + opt$pseudocount )
  b <- log10 ( b + opt$pseudocount )
  
}

output = ifelse(opt$output == "stdout", "", opt$output)


#********
# BEGIN *
#********

# m <- ( a+b )/2

m <- list()

for (my.col in 1:ncol(a)) {

  tmp.df <- cbind(a[, my.col], b[,my.col])
  m[[my.col]] <- rowMeans(tmp.df[,1:2], na.rm=opt$remove_NAs)

}

m <- data.frame(m)
rownames(m) <- rownames(a)
colnames(m) <- colnames(a)


if ( !is.null(opt$round) ){
  
  m <- round(m, opt$round)
  
}




#*********
# OUTPUT *
#*********

write.table(m, file=output, quote=FALSE, sep="\t")



#*******
# EXIT *
#*******

quit(save="no")

