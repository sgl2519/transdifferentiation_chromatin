#************
# LIBRARIES *
#************

library(ggplot2)


#********
# BEGIN *
#********


# 1. source Rscripts for promoters, pEnhancers and dEnhancers
source("fig_4f.R")
source("fig_4g.R")

x2.promoters$region <- "promoters"
x2.proxEnhancers$region <- "proximal enhancers"
x2.distEnhancers$region <- "distal enhancers"



# 2. merge the three types of regions in a unique df
x2 <- rbind(x2.promoters, x2.proxEnhancers, x2.distEnhancers)
x2$combo <- paste(x2$mark, x2$region, sep = "_")


# 3. read list of activated genes unmarked by H3K27me3
act.genes <- read.table('../expression/257.notExpressed.0h.txt')

# Load H3K27me3 marking
h3k27me3 <- read.table('H3K27me3/QN.merged/ant.del.analysis/H3K27me3.peaks.dynamics.tsv')

h3k27me3$V2 <- gsub(';.*', '', h3k27me3$V2)
h3k27me3$V1 <- gsub('\\..*', '', h3k27me3$V1)

act.genes.premarked <- act.genes[act.genes$V1 %in% h3k27me3[h3k27me3$V2 == "H000" |
                                                              h3k27me3$V2 == "H003", ]$V1, ]
act.genes.unmarked <- act.genes[act.genes$V1 %in% setdiff(act.genes$V1, act.genes.premarked), ]


# 4. keep only activated genes unmarked by H3K27me3
x2 <- x2[x2$gene_id %in% act.genes.unmarked, ]
x2 <- x2[complete.cases(x2), ]
x2$region <- factor(x2$region, levels = c("promoters", "proximal enhancers", "distal enhancers"))
x2 <- x2[x2$mark %in% c("H3K4me1",
                        "H3K4me2",
                        "H3K27ac",
                        "H3K9ac",
                        "H3K4me3",
                        "H3K36me3",
                        "H4K20me1",
                        "H3K9me3",
                        "H3K27me3"), ]
x2$mark <- factor(x2$mark,
                  levels = rev(c("H3K4me1", "H3K4me2",
                                 "H3K27ac", "H3K9ac",
                                 "H3K4me3",
                                 "H3K36me3",
                                 "H4K20me1",
                                 "H3K9me3",
                                 "H3K27me3")))

# 5. define color palette
palette <- c("H3K9ac" = "#af4d85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#d199b9",
             "H3K27me3" = "#1d2976",
             "H3K9me3" = "#a7add4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#4c7027",
             "H3K4me1" = "#e5ab00",
             "H3K4me2" = "#a67c00")


# 6. make plot
pdf("fig_S11a_left.pdf", 
    height = 4, width = 6.5, useDingbats = F)
ggplot(data = x2) +
  geom_violin(aes(x=mark, y=value, fill=mark),
              alpha=.5, width = 1, color = "white") +
  stat_summary(mapping = aes(x = mark, y = value),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median, size = .2) +
  # geom_boxplot(data = x2.unmarked.H3K27me3, aes(x=mark, y=value, fill=mark),
  #              width=0.2, alpha = .8, outlier.shape = NA) +
  # geom_text(data = text.df, aes(y=pos, x=Var1, label = paste0(Freq, "%"))) +
  facet_grid(.~region) +
  guides(fill='none') +
  theme_bw() +
  coord_flip() +
  theme(axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size=10, angle = 45, vjust = .95, hjust = .95),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 12),
        panel.border = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        plot.title = element_text(size = 12, hjust = .5)) +
  scale_fill_manual(values = palette) +
  geom_hline(yintercept = 0, color = "#525252", linetype="dashed") +
  ylab("time of histone marking deposition") +
  scale_y_continuous(labels = c("100 h\nbefore", "gene\nactivation", "100 h\nafter", "never\nmarked"),
                     breaks = c(-100, 0, 100, 170))
