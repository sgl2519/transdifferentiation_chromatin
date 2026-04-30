setwd("../analysis/cCREs/precedence/")

library(reshape2)
library(dplyr)
library(pheatmap)
library(ComplexHeatmap)

marks <- c("H3K4me1", "H3K4me2", "H3K27ac", 
           "H3K9ac", "H3K4me3", "H3K36me3", 
           "H4K20me1", "H3K9me3", "H3K27me3")

hours <- c("H000", "H003", "H006", 
           "H009","H012", "H018", 
           "H024", "H036", "H048", 
           "H072", "H120", "H168")

# Load identifier of cCREs with ATAC-seq support
c <- read.table('cCREs.ATACseq_overlap.tsv')

m <- c[c$V8 != '.', ]
m <- unique(m[, 6:7])

c <- subset(c, !(V7 %in% m$V7))

# Build matrix of combinations
m <- read.table('GRCh38-cCREs.V4_wID.bed')
m <- m[m$V7 %in% c$V7, ]
m <- m[6:7]
colnames(m) <- c('class', 'ID')
m <- data.frame(ID = paste0(rep(m$ID, each = 12), '_', hours),
                class = rep(m$class, each = 12))

for ( i in 1:9 ) {
  tmp <- read.table(paste0(marks[i],
                           '/',
                           marks[i],
                           '.peaks.dynamics.binary.tsv'))
  tmp <- tmp[tmp$V1 %in% c$V7, ]
  colnames(tmp) <- c('cCRE_id', hours)
  
  tmp <- reshape2::melt(tmp)
  
  tmp$ID <- paste0(tmp$cCRE_id, '_', tmp$variable)
  tmp <- tmp[, c(4, 3)]
  colnames(tmp) <- c('ID', marks[i])
  m <- merge(m, tmp, by = 'ID', all.x = T)
  
}

rm(tmp)

m[is.na(m)] <- 0

m$superclass <- 'CA/TF'
m[m$class == 'PLS',]$superclass <- 'PLS'
m[m$class == 'pELS',]$superclass <- 'pELS'
m[m$class == 'dELS',]$superclass <- 'dELS'

m$class <- NULL

# Represent cumulative proportion of combinations -- discarding gene body marks
t <- m[, c(2:6, 9:10)] %>%
  group_by(across(everything())) %>%
  count()

t <- as.data.frame(t)

t <- t[order(t$n, decreasing = T), ]

# Represent absolute number of combinations
rownames(t) <- seq(1, nrow(t), 1)

t <- mutate(t, csum = cumsum(n))

cummul <- t$csum/sum(t$n)*100

t <- t[c(which(cummul <= 90), length(which(cummul <= 90)) + 1), ]
t <- as.data.frame(t)
rownames(t) <- 1:nrow(t)

t[, 1:7] <- t(apply(t[, 1:7], 1, function(a) a*(1:7)))

row_h = rowAnnotation('% cCREs' = anno_barplot(cummul[c(which(cummul <= 90), length(which(cummul <= 90)) + 1)],
                                               gp = gpar(fill = 'black')),
                      number = anno_text(paste0(round(cummul[c(which(cummul <= 90), length(which(cummul <= 90)) + 1)], digits = 2), ' %'), 
                                         gp = gpar(fontsize = 12)))


pdf('fig_S1d.left.pdf',
    height = 6,
    width = 5)

print(Heatmap(as.matrix(t[, c(1:7)]), 
              cluster_rows = FALSE,
              cluster_columns = FALSE,
              width = ncol(t)*unit(5, "mm"), 
              height = nrow(t)*unit(5, "mm"),
              show_row_names = TRUE,
              row_names_gp = gpar(fontface = "italic"),
              right_annotation = row_h,
              row_names_side = "left",
              col = c("0" = 'lightgrey',
                      "4" = "#AF4D85",
                      "3" = "#630039",
                      "5" = "#D199B9",
                      "7" = "#1D2976",
                      "6" = "#A7ADD4",
                      "1" = "#E5AB00",
                      "2" = "#A67C00"),
              heatmap_legend_param = list(title = "Marking",
                                          at = c(0:7),
                                          labels = c('absent', 
                                                     'H3K4me1', 'H3K4me2', 
                                                     'H3K27ac', 'H3K9ac', 
                                                     'H3K4me3', 
                                                     'H3K9me3', 'H3K27me3')),
              border_gp = gpar(col = "black"),
              rect_gp = gpar(col = "white", lwd = 2),
              column_title = paste0(format(nrow(m)/12,big.mark=",", trim=TRUE), 
                                    ' cCREs'))
)
dev.off()
