library(pheatmap)
library(fastcluster)
library(dplyr)
library(alluvial)
library(tidyr)
library(ggalluvial)
library(data.table) 
library(ggplot2)
library(mgsub)

matrix <- read.table("../analysis/cCREs/tHMM/HMM.6.gene.matrix.tsv", header = T)

tmp <- matrix
tmp$peakID <- rownames(matrix)

### Generate the annotation
classif <- read.table("/users/project/encode_005982_no_backup/flagship/cCREs/GRCh38-cCREs.V4_wID.bed")
tmp <- merge(x = tmp, y = classif[,6:7], by.x = "peakID", by.y = "V7", all.x = TRUE)
tmp$ATAC <- rep("N", nrow(tmp))
ATAC <- scan("/users/project/encode_005982_no_backup/flagship/cCREs/ATAC.750_overlap.txt",
             character(), quote = "")
tmp$ATAC[tmp$peakID %in% ATAC] <- "Y"       

annot_row <- data.frame(row.names = tmp[, 1],
                        `ENCODEv4` = as.factor(tmp[, 14]),
                        ATAC = as.factor(tmp[, 15]))

prueba <- t(apply(matrix, 1, function(a) mgsub::mgsub(a,
                                                      pattern = c('1', '2', '4', '5', '6', '3'),
                                                      replacement = c('1', '2', '3', '4', '5', '6'))))
             
prueba <- t(apply(prueba, 1, as.numeric))

matrix <- prueba
tmp <- matrix
ser.obj <- seriate(as.matrix(tmp), method = "Mean")
df <- tmp[get_order(ser.obj, 1), ]

pdf('fig_2b.pdf',
    height = 11,
    width = 5.5)

# Plot grouping by cCRE cateory
df <- data.frame()
for ( i in names(encode_colors) ) {
  tmp <- prueba[rownames(annot_row[annot_row$ENCODEv4 == i, ]), ]
  
  ser.obj <- seriate(as.matrix(tmp), method = "Mean")
  df <- rbind(df, tmp[get_order(ser.obj, 1), ])
}

colnames(df) <- c("H000", "H003", "H006",
                  "H009", "H012", "H018",
                  "H024", "H036", "H048",
                  "H072", "H120", "H168")

Heatmap(annot_row[rownames(df), 'ENCODEv4', drop = FALSE],
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = FALSE,
        row_split = factor(annot_row[rownames(df), ]$ENCODEv4,
                           levels = names(encode_colors)),
        row_title_rot = 0,
        width = unit(5, "mm"), 
        col = encode_colors,
        border = TRUE,
        use_raster = TRUE) +
  
  Heatmap(df,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          show_row_names = FALSE,
          width = ncol(df)*unit(5, "mm"), 
          col = c("1" = "#bdbdbd",
                  "2" = "#f7d68f",
                  "6" = "black",
                  "3" = "#E6AB02",
                  "4" = "#e67f30",
                  "5" = "#E7298A"),
          border = TRUE,
          use_raster = TRUE,
          raster_by_magick = FALSE)

dev.off()
