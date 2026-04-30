library(plyr)
library(dplyr)
library(ComplexHeatmap)

setwd("../analysis/cCREs/precedence/")

marks <- c("H3K4me1", "H3K4me2",
           "H3K27ac", "H3K9ac", 
           "H3K4me3", 
           "H3K9me3", "H3K27me3")

# Load the matrix of mark appearance
m <- read.table('mark.gain_tp.tsv',
                header = T)
m <- m[, c("cCRE_id",
           "category",
           marks)]

# Load ATAC-seq support
c <- read.table('cCREs.ATACseq_overlap.tsv')
m <- m[m$cCRE_id %in% c[c$V8 != '.', ]$V7, ]

# Remove cases that do not gain marking
m <- m[apply(m[, 3:(length(marks) + 2)] != 0, 1, any), ]

# Homogenize the order between cCREs
m[, 3:(length(marks) + 2)] <- t(apply(m[, 3:(length(marks) + 2)], 1, function(a) {
  a[a==0] <- NA
  dplyr::dense_rank(a) }))

tmp <- m[, 2:(length(marks) + 2)] %>%
  group_by_all() %>%
  dplyr::summarise(count = n()) %>%
  arrange(-count) %>%
  as.data.frame()

tmp <- ddply(tmp,
             .(category),
             transform,
             prop = count/sum(count)*100)

tmp1 <- tmp[, 2:(length(marks) + 2)] %>%
  group_by(dplyr::across(all_of(marks))) %>%
  dplyr::summarise(total_count = sum(count)) %>%
  arrange(-total_count) %>%
  as.data.frame()

tmp1$total_proportion <- tmp1$total_count/sum(tmp1$total_count)*100

# Add cumulative sum
cummul <- cumsum(tmp1$total_proportion)

pdf('fig_S2a.pdf',
    height = 9,
    width = 4)

row_h = rowAnnotation('% cCREs' = anno_barplot(cummul[c(which(cummul <= 75), length(which(cummul <= 75)) + 1)],
                                               gp = gpar(fill = 'black')),
                      number = anno_text(paste0(round(cummul[c(which(cummul <= 75), length(which(cummul <= 75)) + 1)], digits = 2), ' %'), 
                                         gp = gpar(fontsize = 12)))

plot <- apply(tmp1[c(which(cummul <= 75), length(which(cummul <= 75)) + 1), 1:(length(marks))], 2, as.character)
rownames(plot) <- 1:nrow(plot)
Heatmap(plot,
        column_title = paste0('ATAC+: ', 
                              format(sum(tmp1$total_count),
                                     big.mark=",", trim=TRUE), 
                              ' cCREs'),
        width = ncol(tmp1[c(which(cummul <= 75), length(which(cummul <= 75)) + 1), ])*unit(5, "mm"), 
        height = nrow(tmp1[c(which(cummul <= 75), length(which(cummul <= 75)) + 1), ])*unit(5, "mm"),
        right_annotation = row_h,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_row_names = TRUE,
        row_names_side = "left",
        row_names_gp = gpar(fontface = "italic"),
        col = c("1" = "#49006A",
                "2" = "#AE017E",
                "3" = "#F768A1",
                "4" = "#FCC5C0"),
        na_col = "lightgrey",
        border = TRUE,
        rect_gp = gpar(col = "white", 
                       lwd = 1),
        heatmap_legend_param = list(title = "Deposition\norder", 
                                    at = c("1", "2", "3", "4"),
                                    labels = c("1st", "2nd", "3rd", "4th"))
)

dev.off()
