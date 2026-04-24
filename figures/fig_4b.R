#************
# LIBRARIES *
#************

library(ggplot2)
library(cowplot)
library(scales)
library(circlize)
library(ComplexHeatmap)
library(reshape2)
library(dplyr)

#************
# FUNCTIONS *
#************

function1 <- function(gene) {
  
  # 1. collect expression and marks' 
  # time-series profiles for the selected gene
  df <- data.frame(stringsAsFactors = F)
  
  ## 1.1. expression
  tmp <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
  tmp <- tmp[gene, ]
  tmp$type <- "expression"
  # tmp[1, ] <- (2^(tmp[1, ]) -1) # get raw TPMs
  df <- rbind(df, tmp)
  
  ## 1.2. marks
  for ( i in 1:9 ) {
    
    tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".matrix.tsv"))
    tmp <- tmp[gene, ]
    tmp$type <- marks[i]
    df <- rbind(df, tmp)
    
  }
  
  # 2. rename colnames
  colnames(df)[1:12] <- hours
  
  return(df)
  
}

function2 <- function(gene) {
  
  # 1. collect expression and marks' 
  # time-series profiles for the selected gene
  df <- data.frame(stringsAsFactors = F)
  
  ## 1.1. expression
  tmp <- read.table("expression/QN.merged/selected.genes.rep.2.3.after.QN.merged.tsv", h=T, sep="\t")
  tmp <- tmp[gene, ]
  tmp$type <- "expression"
  # tmp[1, ] <- (2^(tmp[1, ]) -1) # get raw TPMs
  df <- rbind(df, tmp)
  
  ## 1.2. marks
  for ( i in 1:9 ) {
    
    tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".matrix.after.QN.merged.tsv"))
    tmp <- tmp[gene, ]
    tmp$type <- marks[i]
    df <- rbind(df, tmp)
    
  }
  
  # 2. rename colnames
  colnames(df)[1:12] <- hours
  
  return(df)
  
}

#********
# BEGIN *
#********

# 1. set working directory
setwd("../analysis/all.marks")

# 2. the marks we're analyzing
marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3", "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3")

# 3. define color palette
palette <- c("H3K9ac" = "#af4d85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#d199b9",
             "H3K27me3" = "#1d2976",
             "H3K9me3" = "#a7add4",
             "H3K36me3" = "#7fbc41",
             "H4K20me1" = "#4c7027",
             "H3K4me1" = "#e5ab00",
             "H3K4me2" = "#a67c00",
             "expression" = "black")

# 4. vector of hours
hours <- c("H000", "H003", "H006", "H009", "H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")

# 5. define order of marks and expression
my.order <- c("H3K4me1", "H3K4me2", "H3K27ac",
              "expression", "H3K9ac", "H3K4me3",
              "H3K36me3", "H4K20me1", "H3K9me3", "H3K27me3")

#--------------------------------
# the selected gene - CCL2
gene = "ENSG00000108691"

# 4. collect expression and marks' 
# time-series profiles for the selected gene
x <- data.frame(stringsAsFactors = F)

## 4.1. expression
tmp <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
tmp <- tmp[gene, ]
tmp$type <- "expression"
x <- rbind(x, tmp)

## 4.2. marks
for ( i in 1:7 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".matrix.tsv"))
  tmp <- tmp[gene, ]
  tmp$type <- marks[i]
  x <- rbind(x, tmp)
  
}

# 5. rename colnames
hours <- c("H000", "H003", "H006", "H009", "H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")
colnames(x)[1:12] <- hours

# 6. prepare df for plot
x.melt <- melt(x)
x.melt$type <- factor(x.melt$type, levels = c("H3K4me1", "H3K4me2", "H3K27ac",
                                              "expression", "H3K9ac", "H3K4me3",
                                              "H3K36me3", "H4K20me1"))

x.melt$type2 <- ifelse(x.melt$type == "expression", "expression", "mark")
f2 <- max(x.melt[x.melt$type=="expression", "value"]) / max(x.melt[x.melt$type!="expression", "value"])
x.melt[x.melt$type2 == "mark", "value"] <- x.melt[x.melt$type2 == "mark", "value"] * f2

# 13. lineplot
pdf("fig_4b_CCL2.pdf.pdf",
    width = 4.5, height = 3.5)
ggplot(x.melt, aes(x=variable, y=value, group = type, color = type)) +
  geom_line(data=x.melt %>% filter(type2=="expression"), size = 1.5) +
  geom_line(data=x.melt %>% filter(type2=="mark")) +
  scale_color_manual(values = palette) +
  theme_bw() +
  labs(title = gene) +
  guides(color = F) +
  labs ( title = "CCL2") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15, angle = 90, vjust = .5, hjust = .95),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=15, hjust = .5)) +
  ylab("expression (TPM)") +
  scale_y_continuous(
    limits = c(0, 11.98406), #set max break for expression log2(4050+1)    
    breaks = c(0, 2.584963, 5.672425, 8.968667, 11.98406),
    labels = c("0", "5", "50", "500", "4,050"),
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./f2, name = "histone marks (signal)"))
dev.off()

#--------------------------------
# the selected gene - ILK
gene = "ENSG00000166333"

# 4. collect expression and marks' 
# time-series profiles for the selected gene
x <- data.frame(stringsAsFactors = F)

## 4.1. expression
tmp <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
tmp <- tmp[gene, ]
tmp$type <- "expression"
x <- rbind(x, tmp)

## 4.2. marks
for ( i in 1:7 ) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/", marks[i], ".matrix.tsv"))
  tmp <- tmp[gene, ]
  tmp$type <- marks[i]
  x <- rbind(x, tmp)
  
}

# 5. rename colnames
hours <- c("H000", "H003", "H006", "H009", "H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")
colnames(x)[1:12] <- hours

# 6. prepare df for plot
x.melt <- melt(x)
x.melt$type <- factor(x.melt$type, levels = c("H3K4me1", "H3K4me2", "H3K27ac",
                                              "expression", "H3K9ac", "H3K4me3",
                                              "H3K36me3", "H4K20me1"))

x.melt$type2 <- ifelse(x.melt$type == "expression", "expression", "mark")
f2 <- max(x.melt[x.melt$type=="expression", "value"]) / max(x.melt[x.melt$type!="expression", "value"])
x.melt[x.melt$type2 == "mark", "value"] <- x.melt[x.melt$type2 == "mark", "value"] *  f2 # conversion factor is max expression / max mark

# 13. lineplot
pdf("fig_4b_ILK.pdf",
    width = 4.5, height = 3.5)
ggplot(x.melt, aes(x=variable, y=value, color = type, group=type)) +
  geom_line(data=x.melt %>% filter(type2=="mark")) +
  geom_line(data=x.melt %>% filter(type2=="expression"), size = 1.5) +
  scale_color_manual(values = palette) +
  geom_hline(yintercept = 1,
             linetype = 'dashed',
             color = 'lightgrey') +
  theme_bw() +
  guides(color = F) +
  labs ( title = "ILK") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=15, angle = 90, vjust = .5, hjust = .95),
        panel.border = element_rect(color="black"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=15, hjust = .5)) +
  ylab("expression (TPM)") +
  scale_y_continuous(
    # Add a second axis and specify its features
    breaks = c(log2(16+1), log2(32+1), log2(64+1), log2(128+1), log2(256+1)),
    labels = c("16", "32", "64", "128","256"),
    sec.axis = sec_axis(~./f2, name = "histone marks (signal)"))

dev.off()
