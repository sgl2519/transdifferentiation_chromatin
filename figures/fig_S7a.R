#************
# LIBRARIES *
#************

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(cowplot)
library(reshape2)
library(dplyr)
library(tidyr)
library(ggcorrplot)
library(ggsignif)
library(ggpubr)
library(grid)


#************
# FUNCTIONS *
#************

# compute correlation between expression and mark across time
function1 <- function(mark, x="pearson") {
  
  # 1. read mark matrix
  m <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                  h=T, sep="\t")
  
  # 2. check rownames are the same in expression and mark matrix
  stopifnot(identical(rownames(e), rownames(m)))
  
  # 3. compute correlation for each gene between expression and mark profiles
  y <- diag(cor(t(e), t(m), method=x))
  
  # 4. prepare df
  df <- data.frame(cc=y)
  rownames(df) <- rownames(e)
  df$mark <- mark
  df$type <- "across\ntime"
  df$tp <- NA
  
  return(df)
  
}


# steady state correlations between expression and mark
function2 <- function(mark, x="pearson") {
  
  # 1. read mark matrix
  m <- read.table(paste0(mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                  h=T, sep="\t")
  
  # 2. check rownames are the same in expression and mark matrix
  stopifnot(identical(rownames(e), rownames(m)))
  
  # 3. compute correlation across genes between expression and mark values
  y <- diag(cor(e, m, method=x))
  
  # 4. prepare df
  df <- data.frame(cc=y)
  df$mark <- mark
  df$type <- "steady\nstate"
  df$tp <- hours[rownames(df)]
  
  return(df)
  
}



#********
# BEGIN *
#********


# 1. set working directory
setwd("../analysis/all.marks/")

# 2. read expression matrix
e <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv",
                h=T, sep="\t")

# 3. sort rownames of expression matrix
new.order <- sort(rownames(e))
e <- e[new.order, ]


# 4. the marks we're analyzing
marks <- c("H3K27ac", "H3K9ac", "H4K20me1","H3K36me3", "H3K4me3", "H3K4me1",
           "H3K4me2", "H3K9me3", "H3K27me3")

# 5. compute cc for all marks across time-points
w <- data.frame()

for (mark in marks) {
  
  w <- rbind(w, function1(mark=mark))
  
}

w$mark <- factor(w$mark, levels = marks)

w.mean <- w %>%
  group_by(mark) %>%
  summarise(mean = mean(cc, na.rm = T), median = median(cc, na.rm = T))


# 6. compute steady state cc for all marks
hours <- c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168)
names(hours) <- c("H000", "H003", "H006", "H009", "H012",
                  "H018", "H024", "H036", "H048", "H072",
                  "H120", "H168")


z <- data.frame()
for (mark in marks) {
  
  z <- rbind(z, function2(mark=mark))
  
}

z$mark <- factor(z$mark, levels = marks)

z.mean <- z %>%
  group_by(mark) %>%
  summarise(mean = mean(cc, na.rm = T), median = median(cc, na.rm = T))
z$coord <- rep(1:12, 9)


# 7. define color palette
palette <- c("H3K9ac" = "#af4d85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#d199b9",
             "H3K27me3" = "#1d2976",
             "H3K9me3" = "#a7add4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#4c7027",
             "H3K4me1" = "#e5ab00",
             "H3K4me2" = "#a67c00")


# 8. make plot
pdf("../figures/fig_S7a.pdf",
    width = 16, height = 4.5, useDingbats = F)
ggplot() +
  geom_violin(data = w, 
              aes(x=6.5, y=cc, fill=mark), 
              alpha=.3, colour="white", width = 8) +
  geom_boxplot(data = w, 
               aes(x=6.5, y=cc, fill=mark),
               alpha = .4, outlier.shape = NA, width = 3) +
  geom_point(data = z,
             aes(x=coord, y=cc, size = tp, fill = mark),
             shape = 21) +
  guides(fill=F, size = guide_legend(title = "hours", nrow = 1)) +
  ylab(expression(paste("Pearson ", italic("r")))) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=23),
        axis.text.y = element_text(size=23),
        axis.text.x = element_blank(),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 23),
        legend.title = element_text(size = 23),
        strip.background.x = element_blank(),
        strip.text.x = element_text(size=23),
        panel.spacing = unit(0.1, "lines")) +
  scale_fill_manual(values = palette) +
  facet_grid(~mark) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#737373") +
  labs(size = "hours") +
  ylim(-1, 1) +
  scale_size_continuous(breaks = c(0, 3, 6, 9, 12, 18, 24, 36, 48, 72, 120, 168))

dev.off()


