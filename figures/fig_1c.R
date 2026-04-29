setwd('../analysis/cCREs/precedence/')

library(ggplot2)
library(reshape2)
library(data.table)
library(dplyr)

marks <- c("H3K4me1", "H3K4me2", "H3K27ac",
           "H3K9ac", "H3K36me3", "H4K20me1", 
           "H3K4me3", "H3K9me3", "H3K27me3")

hours <-  c("H000", "H003", "H006", 
            "H009","H012", "H018", 
            "H024", "H036", "H048", 
            "H072", "H120", "H168")

x <- read.table('GRCh38-cCREs.V4_wID.bed')
x <- x[c(6, 7)]

# Load ATAC support
a <- as.data.frame(fread('cCREs.ATACseq_overlap.tsv'))

## Extract cCREs with ATAC-support
a <- a[a$V8 != '.', ]

# Load presence of peaks per time point
n <- read.table(paste0(marks[1],
                       "/",
                       marks[1],
                       ".peaks.dynamics.binary.tsv.bz2"), stringsAsFactors = F)
m <- merge(x, n,
           by.x = c('V7'),
           by.y = 'V1',
           all.x = T)
m[is.na(m)] <- 0
colnames(m) <- c('cCRE_id', 'class', hours)
m <- melt(m)
m$ATAC <- 0
m[m$cCRE_id %in% a$V7, ]$ATAC <- 1
m <- m[, c(2:5)]

m <- m %>% group_by_all() %>% summarise(count = n()) %>% as.data.frame()

colnames(m) <- c('class', 'time', 'presence', 'ATAC', 'count')
m$mark <- marks[1]

for ( i in 2:9 ) {
  n <- read.table(paste0(marks[i],
                         "/",
                         marks[i],
                         ".peaks.dynamics.binary.tsv.bz2"), stringsAsFactors = F)
  tmp <- merge(x, n,
               by.x = c('V7'),
               by.y = 'V1',
               all.x = T)
  tmp[is.na(tmp)] <- 0
  colnames(tmp) <- c('cCRE_id', 'class', hours)
  tmp <- melt(tmp)
  tmp$ATAC <- 0
  tmp[tmp$cCRE_id %in% a$V7, ]$ATAC <- 1
  tmp <- tmp[, c(2:5)]
  
  tmp <- tmp %>% group_by_all() %>% summarise(count = n()) %>% as.data.frame()
  
  colnames(tmp) <- c('class', 'time', 'presence', 'ATAC', 'count')
  tmp$mark <- marks[i]
  
  m <- rbind(m, tmp)
}

# Add total number of cCREs
num <- as.data.frame(table(x$V6))
m <- merge(m, num, by.x = 'class', by.y = 'Var1', all.x = T)
m$proportion <- m$count/m$Freq

# Add total number of cCREs dividing by ATAC support/not-support
a_num <- a[, 6:7] %>% unique() %>% group_by(V6) %>% summarise(count = n()) %>% as.data.frame() %>% cbind(data.frame(ATAC = 1))
a_num <- rbind(a_num,
               data.frame(V6 = a_num$V6,
                          count = num$Freq - a_num$count,
                          ATAC = 0))
m <- merge(m, a_num, by.x = c('class', 'ATAC'), by.y = c('V6', 'ATAC'), all.x = T)
m$proportion_ss <- m$count.x/m$count.y


m$mark <- factor(m$mark, levels = marks)

m$class <- factor(m$class, levels = c('PLS',
                                      'pELS',
                                      'dELS',
                                      'CA-H3K4me3',
                                      'CA',
                                      'CA-CTCF',
                                      'CA-TF',
                                      'TF'))

# Remove gene body marks
pdf('fig_1c.pdf',
    height = 4,
    width = 15)

# Plot all CA/TFs aggregated
plot <- m
plot$new_class <- 'CA/TF'
plot[plot$class == 'PLS', ]$new_class <- 'PLS'
plot[plot$class == 'pELS', ]$new_class <- 'pELS'
plot[plot$class == 'dELS', ]$new_class <- 'dELS'

plot <- plot %>%
  group_by(new_class, ATAC, time, presence, mark) %>%
  dplyr::summarise(new_count.x = sum(count.x),
                   new_count.y = sum(count.y)) %>%
  mutate(new_proportion_ss = new_count.x/new_count.y) %>%
  as.data.frame()

plot$mark <- factor(plot$mark, levels = marks)

plot$new_class <- factor(plot$new_class, levels = c('PLS',
                                                    'pELS',
                                                    'dELS',
                                                    'CA/TF'))

ggplot(data = plot[(plot$presence == 1) & (plot$mark %in% marks[c(1:4, 7:9)]), ],
       aes(x = time,
           y = new_proportion_ss,
           color = mark)) +
  geom_line(aes(group = mark)) +
  geom_point() +
  facet_grid(ATAC ~ new_class,
             labeller = labeller(ATAC = c('0' = paste0(format(sum(unique(plot[plot$ATAC == 0, ]$new_count.y)),
                                                              big.mark = ",",
                                                              scientific = FALSE), 
                                                       " cCREs\nnot supported\nby ATAC-seq"),
                                          '1' = paste0(format(sum(unique(plot[plot$ATAC == 1, ]$new_count.y)),
                                                              big.mark = ",",
                                                              scientific = FALSE), 
                                                       " cCREs\n supported\nby ATAC-seq")))) +
  force_panelsizes(rows = unit(rep(2.5, 2) + 0.1, "cm"), 
                   cols = unit(rep(3.75, 4), "cm")) +
  scale_color_manual(values = c("H3K9ac" = "#AF4D85",
                                "H3K27ac" = "#630039",
                                "H3K4me3" = "#D199B9",
                                "H3K27me3" = "#1D2976",
                                "H3K9me3" = "#A7ADD4",
                                "H3K36me3" = "#7FBC41",
                                "H4K20me1" = "#4C7027",
                                "H3K4me1" = "#E5AB00",
                                "H3K4me2" = "#A67C00")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     name = "Proportion of cCREs overlapping\neach hPTM per category",
                     sec.axis = sec_axis(~ ., name=""),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab('') +
  geom_text(data = unique(plot[, c('new_class', 'ATAC', 'new_count.y')]),
            aes(label = paste0(gsub(' ', '',
                                    format(new_count.y,
                                           big.mark = ",",
                                           scientific = FALSE)), ' ', new_class, 's')),
            x = 1,
            y = 1.15,
            color = 'black',
            size = 3.5,
            hjust = 0) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1.20)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(), 
        axis.line.y.right = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.position = 'bottom')

ggplot(data = plot[(plot$presence == 1) & (plot$mark %in% marks[c(1:4, 7, 9)]), ],
       aes(x = time,
           y = new_proportion_ss,
           color = mark)) +
  geom_line(aes(group = mark)) +
  geom_point() +
  facet_grid(ATAC ~ new_class,
             labeller = labeller(ATAC = c('0' = paste0(format(sum(unique(plot[plot$ATAC == 0, ]$new_count.y)),
                                                              big.mark = ",",
                                                              scientific = FALSE), 
                                                       "\nATAC- cCREs"),
                                          '1' = paste0(format(sum(unique(plot[plot$ATAC == 1, ]$new_count.y)),
                                                              big.mark = ",",
                                                              scientific = FALSE), 
                                                       "\nATAC+ cCREs")))) +
  force_panelsizes(rows = unit(rep(2.5, 2) + 0.1, "cm"), 
                   cols = unit(rep(3.75, 4), "cm")) +
  scale_color_manual(values = c("H3K9ac" = "#AF4D85",
                                "H3K27ac" = "#630039",
                                "H3K4me3" = "#D199B9",
                                "H3K27me3" = "#1D2976",
                                "H3K9me3" = "#A7ADD4",
                                "H3K36me3" = "#7FBC41",
                                "H4K20me1" = "#4C7027",
                                "H3K4me1" = "#E5AB00",
                                "H3K4me2" = "#A67C00")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     name = "Proportion of cCREs overlapping\neach hPTM per category",
                     sec.axis = sec_axis(~ ., name=""),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab('') +
  geom_text(data = unique(plot[, c('new_class', 'ATAC', 'new_count.y')]),
            aes(label = paste0(gsub(' ', '',
                                    format(new_count.y,
                                           big.mark = ",",
                                           scientific = FALSE)), ' ', new_class, 's')),
            x = 1,
            y = 1.15,
            color = 'black',
            size = 3.5,
            hjust = 0) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 1.20)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank(), 
        axis.line.y.right = element_blank(),
        strip.background =element_rect(fill="white"),
        legend.position = 'bottom')

dev.off()
