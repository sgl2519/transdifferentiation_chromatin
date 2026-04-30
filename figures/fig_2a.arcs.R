library(devtools)
library(arcdiagram)
library(igraph)
library(dplyr)
library(ggplot2)
library(reshape2)
library(cowplot)

## Generate function to build the network
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

rescale <- function(x){
  return((x -min(x))/(max(x) - min(x)))
}

# Load the matrix containing the states over time for each one of our analysed cCREs
matrix <- read.table("../analysis/cCREs/tHMM/HMM.6.gene.matrix.tsv", header = T)
matrix$peakID <- rownames(matrix)

ids <- read.table("GRCh38-cCREs.V4_wID.bed", header = F)

matrix <- merge(x = matrix, y = ids[,6:7], by.x = "peakID", by.y = "V7", all.x = TRUE)

# Calculate transitions
tmp <- matrix

# Obtain the transition matrix for this subset of cCREs
transition_m <- table(c(as.matrix(tmp[,2:13])[,-ncol(as.matrix(tmp[,2:13]))]), c(as.matrix(tmp[,2:13])[,-1]))
transition_m <- transition_m / rowSums(transition_m)
transition_m <- as.data.frame.matrix(transition_m)

# Load state names
rownames(transition_m) <- c("a: Absent",
                            "b: Basal",
                            "f: Bivalent",
                            "c: Intermediate",
                            "d: Prominent",
                            "e: Strong")

colnames(transition_m) <- c("a: Absent",
                            "b: Basal",
                            "f: Bivalent",
                            "c: Intermediate",
                            "d: Prominent",
                            "e: Strong")

# Generate transition matrix
transition_m <- melt(as.matrix(transition_m), value.name = "p", varnames=c('start', 'end'))
transition_m$p <- round(transition_m$p, 2)

## Remove transitions between the same state
transition_m <- transition_m[transition_m$start != transition_m$end,]
transition_m <- transition_m[transition_m$p != 0,]

## Generate the palette
p <- ggplot(transition_m, aes(x=paste(start, end), y=p, fill=p )) +
  geom_bar(stat="identity") +
  scale_fill_distiller(palette = "Reds", direction = 1, 
                       limits = c(0, max(transition_m$p))) +
  labs(fill = "Transition probability") +
  theme(legend.text = element_text(size=18),
        legend.title = element_text(size=18))

g <- ggplot_build(p)
palette.df <- as.data.frame(as.character(g$data[[1]]["fill"]$fill))

transition_m <- cbind(transition_m, palette = palette.df)

colnames(transition_m) <- c("start", "end", "p", "palette")
## Add number of states suffering each transition
#transition_m$n <- round(nrow(tmp)*11*transition_m$p)

# Sort the states so that the one with higher transition probability are close together

transition_m$start  <- factor(transition_m$start , levels = c("a: Absent",
                                                              "b: Basal",
                                                              "c: Intermediate",
                                                              "d: Prominent",
                                                              "e: Strong",
                                                              "f: Bivalent"))

transition_m$type <- ifelse(as.character(transition_m$start) < as.character(transition_m$end), "increasing", "decreasing")

## Plot the diagram
pdf(paste0("fig_2a.arcs.increasing.pdf"),
    height = 7, width = 8.5)
arcDiagram(as.matrix(transition_m[transition_m$type == "increasing", 1:2]),
           sorted = T,
           lwd = (as.numeric(transition_m[transition_m$type == "increasing",]$p)+0.002)*300,
           #lwd = as.numeric(transition_m[transition_m$type == "increasing",]$n/1100),
           col = transition_m[transition_m$type == "increasing",]$palette,
           mar = c(10.5,0,2,0))
dev.off()


pdf(paste0("fig_2a.arcs.decreasing.pdf"),
    height = 7, width = 8.5)
arcDiagram(as.matrix(transition_m[transition_m$type == "decreasing", 1:2]),
           sorted = T,
           lwd = (as.numeric(transition_m[transition_m$type == "decreasing",]$p)+0.002)*300,
           #lwd = as.numeric(transition_m[transition_m$type == "decreasing",]$n/1100),
           col = transition_m[transition_m$type == "decreasing",]$palette,
           mar = c(10.5,0,2,0))
dev.off()
