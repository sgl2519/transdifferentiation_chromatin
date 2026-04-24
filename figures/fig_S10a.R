# this script plots correlations between expression
# and the nine histone marks for the two sets 
# of 257 activated and 629 already expressed genes
# we plot correlations for all 
# promoters/gene bodies
# proximal regions upstream first TSS
# distal enhancers (ABC blood genewise closest)
# coordinates shifts of distal enhancers (5/10/20/50/100 Kb & 1 Mb)



#************
# LIBRARIES *
#************

library(ggplot2)
library(reshape2)
library(cowplot)



#********
# BEGIN *
#********

# 1. set working directory
setwd("../analysis/")

marks <- c("H3K4me1", "H3K4me2", 
           "H3K27ac", "H3K9ac",
           "H3K4me3",
           "H3K36me3", "H4K20me1",
           "H3K27me3", "H3K9me3")

shift <- c("5Kb", "10Kb", "20Kb", 
           "50Kb", "100Kb", "1Mb")

# 2. define color palette
palette <- c("H3K9ac" = "#AF4D85",
             "H3K27ac" = "#630039",
             "H3K4me3" = "#D199B9",
             "H3K9me3" = "#A7ADD4",
             "H3K27me3" = "#1D2976",
             "H3K36me3" = "#7FBC41",
             "H4K20me1" = "#4C7027",
             "H3K4me1" = "#E5AB00",
             "H3K4me2" = "#A67C00")

act.genes <- read.table('all.marks/expression/257.notExpressed.0h.txt')
alr.expressed.genes <- read.table('all.marks/expression/629.expressed.0h.txt')

# 3. initialize storing objects
m.257 <- data.frame(stringsAsFactors = F)
m.629 <- data.frame(stringsAsFactors = F)

# 4.1. read expression matrix
expression <- read.table('all.marks/expression/QN.merged/expression.matrix.tsv')

# 4.2. for each mark, load the QN matrices at promoters/gene bodies to perform correlations
for (j in c('257 activated genes', 
            '629 already expressed genes')) {
  if ( j == '257 activated genes' ) {
    for ( i in 1:9 ) {
      tmp <- read.table(paste0('all.marks/',
                               marks[i],
                               '/QN.merged/',
                               marks[i],
                               '.matrix.after.QN.merged.tsv'))
      
      # Subset genes of interest
      tmp <- cbind(tmp[act.genes$V1, ], expression[act.genes$V1, ])
      
      # Perform correlation for each gene
      tmp$corr <- apply(tmp, 1, function(a) cor.test(a[1:12], a[13:24],
                                                     method = 'pearson')$estimate)
      
      m.257 <- rbind(m.257,
                     data.frame(corr = tmp$corr,
                                mark = marks[i],
                                region = 'promoters.gb',
                                set = j))
                  
    }
  } else {
    for ( i in 1:9 ) {
      tmp <- read.table(paste0('all.marks/',
                               marks[i],
                               '/QN.merged/',
                               marks[i],
                               '.matrix.after.QN.merged.tsv'))
      
      # Subset genes of interest
      tmp <- cbind(tmp[alr.expressed.genes$V1, ], expression[alr.expressed.genes$V1, ])
      
      # Perform correlation for each gene
      tmp$corr <- apply(tmp, 1, function(a) cor.test(a[1:12], a[13:24],
                                                     method = 'pearson')$estimate)
      
      m.629 <- rbind(m.629,
                     data.frame(corr = tmp$corr,
                                mark = marks[i],
                                region = 'promoters.gb',
                                set = j))
                  
    }
  }
  
}

# 4.3. for each mark, load the QN matrices at proximal enhancers to perform correlations
for (j in c('257 activated genes', 
            '629 already expressed genes')) {
  if ( j == '257 activated genes' ) {
    for ( i in 1:9 ) {
      tmp <- read.table(paste0('enhancers/proximal.enhancers.first.TSS/all.marks/',
                               marks[i],
                               '/QN.merged/',
                               marks[i],
                               '.matrix.after.QN.merged.tsv'))
      
      # Subset genes of interest
      tmp <- cbind(tmp[act.genes$V1, ], expression[act.genes$V1, ])
      
      # Perform correlation for each gene
      tmp$corr <- apply(tmp, 1, function(a) cor.test(a[1:12], a[13:24],
                                                     method = 'pearson')$estimate)
      
      m.257 <- rbind(m.257,
                     data.frame(corr = tmp$corr,
                                mark = marks[i],
                                region = 'proximal.enhancers.first.TSS',
                                set = j))
      
    }
  } else {
    for ( i in 1:9 ) {
      tmp <- read.table(paste0('enhancers/proximal.enhancers.first.TSS/all.marks/',
                               marks[i],
                               '/QN.merged/',
                               marks[i],
                               '.matrix.after.QN.merged.tsv'))
      
      # Subset genes of interest
      tmp <- cbind(tmp[alr.expressed.genes$V1, ], expression[alr.expressed.genes$V1, ])
      
      # Perform correlation for each gene
      tmp$corr <- apply(tmp, 1, function(a) cor.test(a[1:12], a[13:24],
                                                     method = 'pearson')$estimate)
      
      m.629 <- rbind(m.629,
                     data.frame(corr = tmp$corr,
                                mark = marks[i],
                                region = 'proximal.enhancers.first.TSS',
                                set = j))
      
    }
  }
  
}

# 4.4. for each mark, load the QN matrices at distal enhancers to perform correlations
for (j in c('257 activated genes', 
            '629 already expressed genes')) {
  if ( j == '257 activated genes' ) {
    for ( i in 1:9 ) {
      tmp <- read.table(paste0('enhancers/ENCODE-rE2G.blood/all.marks/',
                               marks[i],
                               '/QN.merged/',
                               marks[i],
                               '.matrix.after.QN.merged.tsv'))
      
      # Remove chromosomal region
      rownames(tmp) <- gsub('c.*', '', rownames(tmp))
      
      # Subset genes of interest
      tmp <- cbind(tmp[act.genes$V1, ], expression[act.genes$V1, ])
      tmp <- tmp[complete.cases(tmp), ]
      
      # Perform correlation for each gene
      tmp$corr <- apply(tmp, 1, function(a) cor.test(a[1:12], a[13:24],
                                                     method = 'pearson')$estimate)
      
      m.257 <- rbind(m.257,
                     data.frame(corr = tmp$corr,
                                mark = marks[i],
                                region = 'ENCODE-rE2G.blood.closest',
                                set = j))
      
    }
  } else {
    for ( i in 1:9 ) {
      tmp <- read.table(paste0('enhancers/ENCODE-rE2G.blood/all.marks/',
                               marks[i],
                               '/QN.merged/',
                               marks[i],
                               '.matrix.after.QN.merged.tsv'))
      
      # Remove chromosomal region
      rownames(tmp) <- gsub('c.*', '', rownames(tmp))
      
      # Subset genes of interest
      tmp <- cbind(tmp[alr.expressed.genes$V1, ], expression[alr.expressed.genes$V1, ])
      tmp <- tmp[complete.cases(tmp), ]
      
      # Perform correlation for each gene
      tmp$corr <- apply(tmp, 1, function(a) cor.test(a[1:12], a[13:24],
                                                     method = 'pearson')$estimate)
      
      m.629 <- rbind(m.629,
                     data.frame(corr = tmp$corr,
                                mark = marks[i],
                                region = 'ENCODE-rE2G.blood.closest',
                                set = j))
      
    }
  }
  
}

# 4.5. for each mark, load the QN matrices at shifted distal enhancers to perform correlations
for (j in c('257 activated genes', 
            '629 already expressed genes')) {
  for ( s in shift ) {
    if ( j == '257 activated genes' ) {
      for ( i in 1:9 ) {
        tmp <- read.table(paste0('enhancers/ENCODE-rE2G.blood/shifting/',
                                 s,
                                 '/all.marks/',
                                 marks[i],
                                 '/QN.merged/',
                                 marks[i],
                                 '.matrix.after.QN.merged.tsv'))
        
        # Remove chromosomal region
        rownames(tmp) <- gsub('c.*', '', rownames(tmp))
        
        # Subset genes of interest
        tmp <- cbind(tmp[act.genes$V1, ], expression[act.genes$V1, ])
        tmp <- tmp[complete.cases(tmp), ]
        
        # Perform correlation for each gene
        tmp$corr <- apply(tmp, 1, function(a) cor.test(a[1:12], a[13:24],
                                                       method = 'pearson')$estimate)
        
        m.257 <- rbind(m.257,
                       data.frame(corr = tmp$corr,
                                  mark = marks[i],
                                  region = paste0('ENCODE-rE2G.blood.closest.', s),
                                  set = j))
        
      }
    } else {
      for ( i in 1:9 ) {
        tmp <- read.table(paste0('enhancers/ENCODE-rE2G.blood/shifting/',
                                 s,
                                 '/all.marks/',
                                 marks[i],
                                 '/QN.merged/',
                                 marks[i],
                                 '.matrix.after.QN.merged.tsv'))
        
        # Remove chromosomal region
        rownames(tmp) <- gsub('c.*', '', rownames(tmp))
        
        # Subset genes of interest
        tmp <- cbind(tmp[alr.expressed.genes$V1, ], expression[alr.expressed.genes$V1, ])
        tmp <- tmp[complete.cases(tmp), ]
        
        # Perform correlation for each gene
        tmp$corr <- apply(tmp, 1, function(a) cor.test(a[1:12], a[13:24],
                                                       method = 'pearson')$estimate)
        
        m.629 <- rbind(m.629,
                       data.frame(corr = tmp$corr,
                                  mark = marks[i],
                                  region = paste0('ENCODE-rE2G.blood.closest.', s),
                                  set = j))
        
      }
    }
  }
  
}




sets <- c("promoters.gb", 
          "proximal.enhancers.first.TSS", 
          "ENCODE-rE2G.blood.closest", 
          "ENCODE-rE2G.blood.closest.5Kb", 
          "ENCODE-rE2G.blood.closest.10Kb",
          "ENCODE-rE2G.blood.closest.20Kb",
          "ENCODE-rE2G.blood.closest.50Kb",
          "ENCODE-rE2G.blood.closest.100Kb",
          "ENCODE-rE2G.blood.closest.1Mb")

labels <- c("promoters", "pEnh", "dEnh", "dEnh + 5Kb", "dEnh + 10Kb", "dEnh + 20Kb", "dEnh + 50Kb", "dEnh + 100Kb", "dEnh + 1Mb")

names(labels) <- sets

# 5. plot
lop <- list()

m.257$mark <- factor(m.257$mark, levels = marks)
m.257$region <- factor(m.257$region, levels = sets)

lop[[1]] <- ggplot(data = m.257,
                   aes(x = region, 
                       y = corr, 
                       fill = mark, 
                       alpha = region)) +
  geom_boxplot() +
  facet_wrap(~ mark, nrow = 1) +
  scale_fill_manual(values = palette) +
  guides(fill="none", alpha = 'none') +
  scale_x_discrete(labels = labels) +
  scale_alpha_discrete(range = c(1, 0)) +
  ylab(expression(paste("Pearson ", italic("r")))) +
  labs(title = "257 activated genes") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, hjust=1, vjust=.5, angle=90),
        panel.border = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 15, hjust = .5),
        strip.text = element_text(size = 15),
        legend.position = "bottom",
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")

m.629$mark <- factor(m.629$mark, levels = marks)
m.629$region <- factor(m.629$region, levels = sets)

lop[[2]] <- ggplot(data = m.629,
                   aes(x = region, 
                       y = corr, 
                       fill = mark, 
                       alpha = region)) +
  geom_boxplot() +
  facet_wrap(~ mark, nrow = 1) +
  scale_fill_manual(values = palette) +
  guides(fill="none", alpha = 'none') +
  scale_x_discrete(labels = labels) +
  scale_alpha_discrete(range = c(1, 0)) +
  ylab(expression(paste("Pearson ", italic("r")))) +
  labs(title = "629 already expressed genes") +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.text.x = element_text(size=12, hjust=1, vjust=.5, angle=90),
        panel.border = element_rect(color="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 15, hjust = .5),
        strip.text = element_text(size = 15),
        legend.position = "bottom",
        strip.background = element_blank(),
        panel.spacing = unit(0.1, "lines")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black")

pdf("fig_S10a.pdf", 
    height = 7, width = 15)
plot_grid(plotlist = lop, nrow = 2, align = "v")
dev.off()
