# Author: Ramil Nurtdinov

library(preprocessCore)

file_in <- commandArgs(trailingOnly = TRUE)

if( is.na(file_in)=='TRUE')
{
    print('Bad inline parameters')
    q()
}
        
options(width=150)

data_matrix <- read.table(file=file_in,sep='\t', check.names=FALSE, header=T, row.names=1)
data_matrix <- as.matrix(data_matrix)

new_matrix  <- normalize.quantiles(data_matrix, copy = TRUE)
rownames(new_matrix) <- rownames(data_matrix)
colnames(new_matrix) <- colnames(data_matrix)

write.table(new_matrix, file="", sep='\t', quote=F)

