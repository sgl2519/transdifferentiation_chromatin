#************
# LIBRARIES *
#************

library(dplyr)
library(ggplot2)
library(cowplot)
library(reshape2)
library(ppcor)


#************
# FUNCTIONS *
#************

function1 <- function(mark, set = "257_genes") {
  
  x <- list()
  
  # expression matrix
  x[[1]] <- read.delim("all.marks/expression/QN.merged/expression.matrix.tsv",
                       h=T, sep="\t")
  # promoter matrix
  x[[2]] <- read.delim(paste0("all.marks/", 
                              mark, "/QN.merged/", mark, ".matrix.tsv"),
                       h=T)
  
  # proximal enhancer matrix
  x[[3]] <- read.delim(paste0("enhancers/proximal.enhancers.first.TSS/all.marks/",
                              mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                       h=T)
  
  # distal enhancer matrix
  x[[4]] <- read.delim(paste0("enhancers/ENCODE-rE2G.blood/all.marks/",
                              mark, "/QN.merged/", mark, ".matrix.after.QN.merged.tsv"),
                       h=T)
  rownames(x[[4]]) <- gsub("\\c.*","", rownames(x[[4]]))
  
  
  if (set == "257_genes") {
    
    genes <- read.table("all.marks/expression/257.notExpressed.0h.txt",
                        stringsAsFactors = F) 
    
    
  } else {
    
    genes <- read.table("all.marks/expression/629.expressed.0h.txt",
                        stringsAsFactors = F) 
    
    
  }
  
  colnames(genes) <- "gene_id"
  
  y <- list()
  vars <- c("expression", "promoter", "proxEnh", "distalEnh")
  
  for ( i in 1:4 ){
    
    y[[i]] <- x[[i]][rownames(x[[i]]) %in% genes$gene_id, ]
    y[[i]] <- y[[i]][genes$gene_id, ]
    # y[[i]] <- y[[i]][complete.cases(y[[i]]), ]
    
  }
  
  expression.promoter <- c()
  expression.promoter.a.pr <- c()
  expression.promoter.b.pr <- c()
  expression.promoter.c.pr <- c()
  expression.proxEnh <- c()
  expression.proxEnh.pr <- c()
  expression.distEnh <- c()
  expression.distEnh.pr <- c()
  
  for ( i in 1:nrow(y[[4]]) ) {
    
    y.data <- data.frame(expression = as.numeric(y[[1]][i, ]),
                         promoter = as.numeric(y[[2]][i, ]),
                         proxEnh = as.numeric(y[[3]][i, ]),
                         distEnh = as.numeric(y[[4]][i, ]))
    
    expression.promoter <- c(expression.promoter, 
                             cor(y.data$expression, y.data$promoter, method = "pearson"))
    
    expression.proxEnh <- c(expression.proxEnh,
                            cor(y.data$expression, y.data$proxEnh, method = "pearson"))
    expression.proxEnh.pr <- c(expression.proxEnh.pr, 
                               pcor.test(y.data$expression, y.data$proxEnh, y.data[, c("promoter")])$estimate)
    
    expression.promoter.a.pr <- c(expression.promoter.a.pr, 
                                  pcor.test(y.data$expression, y.data$promoter, y.data[, c("proxEnh")], method = "pearson")$estimate)
    
    
    if (sum(is.na(y.data$distEnh)) != 12) {
      
      expression.distEnh <- c(expression.distEnh, 
                              cor(y.data$expression, y.data$distEnh, method = "pearson"))
      
      expression.promoter.b.pr <- c(expression.promoter.b.pr, 
                                    pcor.test(y.data$expression, y.data$promoter, y.data[, c("distEnh")], method = "pearson")$estimate)
      
      expression.promoter.c.pr <- c(expression.promoter.c.pr, 
                                  pcor.test(y.data$expression, y.data$promoter, y.data[, c("proxEnh", "distEnh")], method = "pearson")$estimate)
      
      expression.distEnh.pr <- c(expression.distEnh.pr,
                                 pcor.test(y.data$expression, y.data$distEnh, y.data[, c("promoter")])$estimate)
      
    }
    
    
  }
  
  plot.m <- data.frame(cc = c(expression.promoter, 
                              expression.promoter.a.pr, 
                              expression.promoter.b.pr,
                              expression.promoter.c.pr,
                              expression.proxEnh,
                              expression.proxEnh.pr,
                              expression.distEnh, 
                              expression.distEnh.pr),
                       group = c(rep("promoter", length(expression.promoter)),
                                 rep("promoter | pEnh", length(expression.promoter.a.pr)),
                                 rep("promoter | dEnh", length(expression.promoter.b.pr)),
                                 rep("promoter | (pEnh, dEnh)", length(expression.promoter.c.pr)),
                                 rep("pEnh", length(expression.proxEnh)),
                                 rep("pEnh | promoter", length(expression.proxEnh.pr)),
                                 rep("dEnh", length(expression.distEnh)),
                                 rep("dEnh | promoter", length(expression.distEnh.pr))))
  
  plot.m$group <- factor(plot.m$group, levels = c("dEnh", "dEnh | promoter",
                                                  "pEnh", "pEnh | promoter",
                                                  "promoter",
                                                  "promoter | dEnh",
                                                  "promoter | pEnh",
                                                  "promoter | (pEnh, dEnh)"))
  
  p <- ggplot(plot.m, aes(x=group, y=cc, fill=group)) +
    geom_violin(alpha=.4, color = "white", scale='width') +
    geom_boxplot(width=.25, alpha = .8, outlier.shape = NA) +
    labs(title = mark) +
    guides(fill='none') +
    ylab(expression(paste("Pearson ", italic("r")))) +
    theme_bw() +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          axis.text.y = element_text(size=13),
          axis.text.x = element_text(size=13, hjust=1, vjust=.5, angle=90),
          panel.border = element_rect(color="black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 15, hjust = .5),
          axis.line = element_line(colour = "black"),
          legend.position = "bottom",
          strip.background = element_blank(),
          strip.text.x = element_text(size = 13, angle = -360),
          panel.spacing = unit(0.1, "lines")) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") +
    scale_fill_manual(values = c("promoter" = "gray",
                                 "promoter | pEnh" = "red",
                                 "promoter | dEnh" = "red",
                                 "promoter | (pEnh, dEnh)" = "red",
                                 "pEnh" = "gray",
                                 "pEnh | promoter" = "red",
                                 "dEnh" = "gray",
                                 "dEnh | promoter" = "red"))
  
  return(p)
  
  
}


#********
# BEGIN *
#********

setwd("../analysis")

marks <- c("H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3")

lop.257 <- list()
lop.629 <- list()

for (l in 1:5) {
  
  lop.257[[l]] <- function1(mark = marks[l])
  lop.629[[l]] <- function1(mark = marks[l], set = "629_genes")
  
  
}

lop.257[2:5] <- lapply(lop.257[2:5], function(x){x <- x + ylab("")})
lop.629[2:5] <- lapply(lop.629[2:5], function(x){x <- x +  ylab("")})

leg.df <- data.frame(group=c("raw cc", "partial cc"),
                     n = c(10, 5))
leg.df$group <- factor(leg.df$group, levels = c("raw cc", "partial cc"))

leg.p <- ggplot(leg.df, aes(x=group, y=n, fill=group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("gray", "red")) +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 12))

leg <- get_legend(leg.p)

pdf("fig_4c.pdf", 
    width = 14, height = 4)
main.p <- plot_grid(plotlist = lop.257, nrow = 1)
plot_grid(main.p, leg, nrow = 2, rel_heights = c(0.9, 0.1))
dev.off()
