#************
# LIBRARIES *
#************
library(arcdiagram)
library(igraph)
library(dplyr)
library(ggplot2)
library(cowplot)


#************
# FUNCTIONS *
#************

# arc diagram
arcDiagram <- function(
  edgelist, sorted=TRUE, decreasing=FALSE, lwd=NULL,
  col=NULL, cex=NULL, col.nodes=NULL, lend=1, ljoin=2, lmitre=1,
  las=2, bg=NULL, mar=c(4,1,3,1))
{
  # ARGUMENTS
  # edgelist:   two-column matrix with edges
  # sorted:     logical to indicate if nodes should be sorted
  # decreasing: logical to indicate type of sorting (used only when sorted=TRUE)
  # lwd:        widths for the arcs (default 1)
  # col:        color for the arcs (default "gray50")
  # cex:        magnification of the nodes labels (default 1)
  # col.nodes:  color of the nodes labels (default "gray50")
  # lend:       the line end style for the arcs (see par)
  # ljoin:      the line join style for the arcs (see par)
  # lmitre:     the line mitre limit fort the arcs (see par)
  # las:        numeric in {0,1,2,3}; the style of axis labels (see par)
  # bg:         background color (default "white")
  # mar:        numeric vector for margins (see par)
  
  # make sure edgelist is a two-col matrix
  if (!is.matrix(edgelist) || ncol(edgelist)!=2)
    stop("argument 'edgelist' must be a two column matrix")
  edges = edgelist
  # how many edges
  ne = nrow(edges)
  # get nodes
  nodes = unique(as.vector(edges))
  nums = seq_along(nodes)
  # how many nodes
  nn = length(nodes)  
  # ennumerate
  if (sorted) {
    nodes = sort(nodes, decreasing=decreasing)
    nums = order(nodes, decreasing=decreasing)
  }
  # check default argument values
  if (is.null(lwd)) lwd = rep(1, ne)
  if (length(lwd) != ne) lwd = rep(lwd, length=ne)
  if (is.null(col)) col = rep("gray50", ne)
  if (length(col) != ne) col = rep(col, length=ne)
  if (is.null(cex)) cex = rep(1, nn)
  if (length(cex) != nn) cex = rep(cex, length=nn)
  if (is.null(col.nodes)) col.nodes = rep("gray50", nn)
  if (length(col.nodes) != nn) col.nodes = rep(col.nodes, length=nn)
  if (is.null(bg)) bg = "white"
  
  # node labels coordinates
  nf = rep(1 / nn, nn)
  # node labels center coordinates
  fin = cumsum(nf)
  ini = c(0, cumsum(nf)[-nn])
  centers = (ini + fin) / 2
  names(centers) = nodes
  
  # arcs coordinates
  # matrix with numeric indices
  e_num = matrix(0, nrow(edges), ncol(edges))
  for (i in 1:nrow(edges))
  {
    e_num[i,1] = centers[which(nodes == edges[i,1])]
    e_num[i,2] = centers[which(nodes == edges[i,2])]
  }
  # max arc radius
  radios = abs(e_num[,1] - e_num[,2]) / 2
  max_radios = which(radios == max(radios))
  max_rad = unique(radios[max_radios] / 2)
  # arc locations
  locs = rowSums(e_num) / 2
  
  # plot
  par(mar = mar, bg = bg)
  plot.new()
  plot.window(xlim=c(-0.025, 1.025), ylim=c(0, 1*max_rad*2))
  # plot connecting arcs
  z = seq(0, pi, l=100)
  for (i in 1:nrow(edges))
  {
    radio = radios[i]
    x = locs[i] + radio * cos(z)
    y = radio * sin(z)
    lines(x, y, col=col[i], lwd=lwd[i], 
          lend=lend, ljoin=ljoin, lmitre=lmitre)
  }
  # add node names
  mtext(nodes, side=1, line=0, at=centers, cex=cex, 
        col=col.nodes, las=las)
}

# count duplicated rows
count.duplicates <- function(DF){
  x <- do.call('paste', c(DF, sep = '\r'))
  ox <- order(x)
  rl <- rle(x[ox])
  cbind(DF[ox[cumsum(rl$lengths)],,drop=FALSE],count = rl$lengths)
  
}


#********
# BEGIN *
#********

# 1. set wd
setwd("../analysis/all.marks")

# 2. read list of 257 not expressed genes
x <- read.table("expression/257.notExpressed.0h.txt",
                stringsAsFactors = F)
colnames(x) <- "gene_id"


# 3. marks 
marks <- c("H3K27ac", "H3K9ac","H3K4me3", "H3K4me1", "H3K4me2", "H3K36me3",
           "H4K20me1", "H3K9me3", "H3K27me3")

# 4. time-points
hours <- c("H000", "H003", "H006", "H009","H012", "H018", "H024", "H036", "H048", "H072", "H120", "H168")


# 5. retrieve time-point of peak appearance of 
# for each histone mark

for ( i in c(1:7, 9)) {
  
  tmp <- read.table(paste0(marks[i], "/QN.merged/ant.del.analysis/", marks[i],
                           ".peaks.dynamics.tsv"), stringsAsFactors = F)
  colnames(tmp) <- c("gene_id", "tp")
  tmp$gene_id <- gsub("\\..*","", tmp$gene_id)
  tmp$tp <- gsub("\\;.*","", tmp$tp)
  colnames(tmp)[2] <- marks[i]
  x <- merge(x, tmp, by = "gene_id", all.x = T)
  
}

# 6. add time-point for expression (aka the 1st time-point in which exp is > 1 TPM)
expression.matrix <- read.table("expression/QN.merged/expression.matrix.tsv", h=T, sep="\t")
expression.matrix <- expression.matrix[rownames(expression.matrix) %in% x$gene_id, ]

m <- data.frame(stringsAsFactors = F)

for (gene in x$gene_id) {
  
  v <- rep(0, 12)
  stn <- as.numeric(expression.matrix[gene, ])
  
  v[which(stn > 1)] <- 1
  
  m <- rbind(m, v)
}

x$expression <- hours[apply(m, 1, function(x){min(which(x>0))})]


# 7. keep genes marked by K27me3 at 0 hours
x <- x[(x$H3K27me3 == "H000" | x$H3K27me3 == "H003") & !(is.na(x$H3K27me3)), ]

# 8. aggregate co-occurring marks
x.melt <- reshape2::melt(x, id.vars = "gene_id")
x.melt <- x.melt[complete.cases(x.melt), ]
x.melt.agg <- aggregate(data=x.melt,variable~gene_id+value,FUN=paste)


# 9. convert lists of co-occurring marks to strings
y <- c()
for (i in 1:nrow(x.melt.agg)) {
  
  stn <- x.melt.agg[i, "variable"]
  stn[[1]] <- factor(stn[[1]], levels = c("H3K4me1", "H3K4me2", "H3K27me3", "H3K4me3",
                                          "H3K27ac", "H3K9ac", "H3K36me3", "H4K20me1", "expression"))
  y <- c(y, paste(sort(stn[[1]]), collapse=";"))
  
}
x.melt.agg$variable <- y



# 11. sort genes by tp, to get full path for every gene across time
x.sorted <- x.melt.agg[with(x.melt.agg, order(gene_id, value)), c("gene_id", "value", "variable")]


# 12. the root groups of marks to which a gene can be assigned 
# are manually defined

root.gs <- c("H3K4me1;H3K4me2;H3K4me3;H3K27me3",
             "H3K9ac;H4K20me1",
             "expression;H3K27ac",
             "H3K36me3")

root.ls <- list()
for ( i in 1:length(root.gs) ){
  root.ls[[i]] <- strsplit(root.gs[i], ";")[[1]]
  
}

# 13. reassign genes to any of the root groups
reassigned.gs <- c()
reassigned.gs <- c()
for ( i in 1:nrow(x.sorted) ){
  if (length(strsplit(x.sorted[i, "variable"], ";")[[1]]) > 1) {
    
    stn <- strsplit(x.sorted[i, "variable"], ";")[[1]]
    
    
    
  } else {
    
    stn <- x.sorted[i, "variable"]
    
  }
  
  
  j_index <- 0
  for ( k in 1:length(root.gs) ){
    
    ints <- (length(intersect(stn, root.ls[[k]])))
    curr_j <- ints/(length(stn) + length(root.ls[[k]]) - ints) 
    if (curr_j > j_index) {
      
      j_index <- curr_j
      assigned_gs <- root.gs[k]
    }
  }
  
  reassigned.gs <- c(reassigned.gs, assigned_gs) 
  
  
}
x.sorted$variable <- reassigned.gs


# 14. reconstruct dynamics path for every gene
# aka order in which expression and marks appear
# co-occurring marks and/or expression are separated by ;
genes <- unique(x.sorted$gene_id)
ss.df <- data.frame(stringsAsFactors = F)

for ( i in 1:length(genes) ) {
  
  vc <- c()
  z <- x.sorted[x.sorted$gene_id == genes[i], ]
  for ( j in 1:nrow(z)) {
    
    vc <- c(vc, z[j, "variable"])
    
  }
  
  vc <- paste(vc, collapse = ">")
  ss.df <- rbind(ss.df, data.frame(gene_id = genes[i], 
                                   ss = vc))
  
}
ss.df$ss <- paste("start", ss.df$ss, sep=">")


# 15. compute pairs of consecutive transitions
v1 <- c()
v2 <- c()

for ( i in 1:nrow(ss.df) ){
  
  stn <- strsplit(ss.df[i, "ss"], ">")[[1]]
  for (j in 1:(length(stn)-1)) {
    
    v1 <- c(v1, stn[j])
    v2 <- c(v2, stn[j+1])
    
  }
  
}

counts.df <- data.frame(v1=v1, v2=v2)
counts.df <- count.duplicates(counts.df)
counts.df$transition <- paste(counts.df$v1, counts.df$v2, sep="_")
counts.df <- counts.df[, c(4,1,2,3)]
colnames(counts.df)[2:4] <- c("i_state", "j_state", "n")
counts.df <- counts.df[counts.df$n >= 10, ] # keep only transitions w/ >= 10 genes 
counts.df$i_state <- as.character(counts.df$i_state)
counts.df$j_state <- as.character(counts.df$j_state)

# 16. reorder combinations of marks
correspondence.v <- c()
names.correspondence.v <- c("start",
                            "H3K4me1;H3K4me2;H3K4me3;H3K27me3",
                            "H3K9ac;H4K20me1",
                            "expression;H3K27ac",
                            "H3K36me3")
correspondence.v <- paste(letters[1:length(names.correspondence.v)], 
                          names.correspondence.v, sep="_")
names(correspondence.v) <- names.correspondence.v

counts.df$i_state <- correspondence.v[counts.df$i_state]
counts.df$j_state <- correspondence.v[counts.df$j_state]


# 17. make plot

pdf("fig_S10d_right.arcs_increase.pdf", height = 5, width = 10)
arcDiagram(as.matrix(counts.df[counts.df$i_state < counts.df$j_state, 2:3]),
           lwd = (as.numeric(counts.df[counts.df$i_state < counts.df$j_stat, "n"]/5)))
dev.off()

pdf("fig_S10d_right.arcs_decrease.pdf", height = 5, width = 10)
arcDiagram(as.matrix(counts.df[counts.df$i_state > counts.df$j_state, 2:3]),
           lwd = (as.numeric(counts.df[counts.df$i_state > counts.df$j_stat, "n"]/5)))
dev.off()
