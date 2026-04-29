library(ggplot2)
library(reshape2)
library(seriation)
library(ComplexHeatmap)
library(tidyr)
library(gridExtra)

marks <- c("H3K4me1", "H3K4me2", "H3K27ac",
           "H3K9ac", "H3K4me3", "H3K36me3", 
           "H4K20me1", "H3K9me3", "H3K27me3")

# Load matrices of peak presence/absence
df <- read.table(paste0('../analysis/cCREs/precedence/',
                        marks[1],
                        '/',
                        marks[1],
                        '.peaks.dynamics.binary_metadata.bz2'),
                 header = TRUE)

df <- df[, -c(2, 3, 4)]
df <- reshape2::melt(df)

for ( i in marks[2:9] ) {
  tmp <- read.table(paste0('../analysis/cCREs/precedence/',
                           i,
                           '/',
                           i,
                           '.peaks.dynamics.binary_metadata.bz2'),
                    header = TRUE)
  
  tmp <- tmp[, -c(2, 3, 4)]
  tmp <- reshape2::melt(tmp)
  
  # Merge with main matrix
  df <- cbind(df, tmp['value'])
  
}

# Add column indicating the number of marks that are concurrently happening
df$num_marks <- apply(df[, 5:13], 1, sum)

# Summarise matrix
plot <- df[, c(2, 3, 14)] %>% 
  group_by(ENCODE4, ATAC.supp, num_marks) %>% 
  dplyr::summarise(sum = n()) %>% 
  as.data.frame()

total_plot <- plot %>% 
  group_by(ATAC.supp, num_marks) %>% 
  dplyr::summarise(sum = sum(sum)) %>% 
  as.data.frame()

plot <- rbind(plot, cbind(data.frame(ENCODE4 = "Total"), total_plot))
plot$ENCODE4 <- factor(plot$ENCODE4, levels = c('Total', 
                                                  'PLS', 'pELS', 'dELS',
                                                  'CA-H3K4me3', 'CA', 'CA-CTCF',
                                                  'CA-TF', 'TF'))
plot$num_marks <- factor(plot$num_marks, levels = rev(sort(unique(plot$num_marks))))

# Find matrix of totals
totals <- plot %>%
  group_by(ENCODE4, ATAC.supp) %>%
  dplyr::summarise(n = sum(sum)) %>%
  as.data.frame()

plot <- merge(plot, totals)

# Prepare plot
colfun <- colorRampPalette(c("#003c30", "#f5f5f5", "#543005"))

# Discard gene body marks
## Add column indicating the number of marks that are concurrently happening
df$num_marks <- apply(df[, c(5:9, 12:13)], 1, sum)

## Summarise matrix
plot <- df[, c(2, 3, 14)] %>% 
  group_by(ENCODE4, ATAC.supp, num_marks) %>% 
  dplyr::summarise(sum = n()) %>% 
  as.data.frame()

total_plot <- plot %>% 
  group_by(ATAC.supp, num_marks) %>% 
  dplyr::summarise(sum = sum(sum)) %>% 
  as.data.frame()

plot <- rbind(plot, cbind(data.frame(ENCODE4 = "Total"), total_plot))
plot$ENCODE4 <- factor(plot$ENCODE4, levels = c('Total', 
                                                'PLS', 'pELS', 'dELS',
                                                'CA-H3K4me3', 'CA', 'CA-CTCF',
                                                'CA-TF', 'TF'))
plot$num_marks <- factor(plot$num_marks, levels = rev(sort(unique(plot$num_marks))))

# Find matrix of totals
totals <- plot %>%
  group_by(ENCODE4, ATAC.supp) %>%
  dplyr::summarise(n = sum(sum)) %>%
  as.data.frame()

plot <- merge(plot, totals)

pdf('fig_1b.pdf',
    height = 3.7,
    width = 4.2)
ggplot(data = plot, 
       aes(x = ENCODE4, 
           y = sum, 
           fill = num_marks)) +
  facet_wrap(ATAC.supp ~ .,
             labeller = labeller(ATAC.supp = c("negative" = paste0(format(sum(plot[plot$ATAC.supp == 'negative', ]$sum/2/12),
                                                                          big.mark = ",",
                                                                          scientific = FALSE), 
                                                                   " ATAC- cCREs"),
                                               "positive" = paste0(format(sum(plot[plot$ATAC.supp == 'positive', ]$sum/2/12),
                                                                          big.mark = ",",
                                                                          scientific = FALSE), 
                                                                   " ATAC+ cCREs")))) +
  geom_bar(stat = "identity",
           position = "fill") +
  scale_y_continuous(labels = scales::percent,
                     breaks = seq(0, 1, by = 0.25)) +
  expand_limits(y = c(0, 1.35),
                x = c(-1.5, 9)) +
  scale_fill_manual(values = c('0' = 'lightgrey', 
                               '1' = "#543005", 
                               '2' = "#8c510a", 
                               '3' = "#bf812d", 
                               '4' = "#aed5d2", 
                               '5' = "#35978f", 
                               '6' = "#01665e", 
                               '7' = "#003c30")) +
  ylab("% cCREs") +
  xlab("") +
  geom_text(data = totals[totals$ENCODE4 != 'Total', ],
            aes(x = ENCODE4, 
                y = 1.01, 
                label = gsub(' ', '',
                             format(n/12,
                               big.mark = ",",
                               scientific = FALSE)),
                fill = NULL),
            angle = 80,
            vjust = 0,
            hjust = 0,
            size = 2.5) +
  guides(fill = guide_legend(title = "Number of\nmarks concurrently\noverlapped")) +
  annotate("text", 
           label = "# cCREs",
           fontface = "bold",
           x = 0, 
           y = 1.1, 
           size = 2.5, 
           colour = "black",
           vjust = 0,
           hjust = 0.5,) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,
                                   face = ifelse(levels(plot$ENCODE4) == "Total","bold","plain")),
        legend.position = "bottom")
dev.off()
