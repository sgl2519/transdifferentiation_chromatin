library(ggplot2)
library(ggsankey)
library(dplyr)
library(stringr)
library(reshape2)

# Set function for finding consecutive elements
find_consecutive_with_values <- function(vec) {
  # Initialize result vectors for counts and values
  counts <- numeric(0)
  values <- numeric(0)
  
  # Start with the first element
  current_value <- vec[1]
  count <- 1
  
  # Loop through the vector starting from the second element
  for (i in 2:length(vec)) {
    if (vec[i] == current_value) {
      # If the same as the previous one, increment the count
      count <- count + 1
    } else {
      # If different, append the count and the value to the results
      counts <- c(counts, count)
      values <- c(values, current_value)
      count <- 1
      current_value <- vec[i]
    }
  }
  
  # Append the last group count and value
  counts <- c(counts, count)
  values <- c(values, current_value)
  
  return(list(counts = counts, values = values))
}


# Load cCREs that gain marking
setwd("../analysis/cCREs/precedence/")

marks <- c("H3K27ac", "H3K9ac", "H3K4me3", 
           "H3K4me1", "H3K4me2", "H3K36me3", 
           "H4K20me1", "H3K9me3", "H3K27me3", "ATAC-seq")

# Load the matrix of mark appearance
m <- read.table('mark.loss_tp.tsv',
                header = T)

# Remove the cases where no active mark appears during transdifferentiation
m <- m[rowSums(m[, c(3:9)]) > 0, ]

# Remove the cases where no active mark appears during transdifferentiation (without considering H3K4me1/2)
#m <- m[rowSums(m[, c(3:5, 8:9)]) > 0, ]
m <- m[rowSums(m[, c(3:5)]) > 0, ]

# Select only those cases where all active marks appear (i.e, are not already present) during the process
#m <- m[rowSums(m[, c(3:5, 8:9)] == 1)==0, , drop = FALSE]
m <- m[rowSums(m[, c(3:5)] == 11)==0, , drop = FALSE]
m <- m[rowSums(m[, c(3:5)] == 12)==0, , drop = FALSE]

# Select cases that when gaining marking, they keep it for a determined (cons_tps) number of time points
cons_tps <- 3

consistent <- list()
for ( i in marks[1:5] ) {
  # Read peak presence and absence matrix
  tmp <- read.table(paste0(i,
                           '/',
                           i,
                           '.peaks.dynamics.binary.tsv'))
  
  # Subset the cCREs from the previous matrix
  tmp <- tmp[tmp$V1 %in% m$cCRE_id, ]
  
  
  # Calculate the number of consecutive time points with marking
  tmp$consecutive <- apply(tmp[, 2:13], 1, function(a) {
    x <- find_consecutive_with_values(as.integer(a))
    if( any(x$counts[x$values == 1] >= cons_tps) ) {
      return('consistent')
    } else {
      return('inconsistent')
    }
    
  })
  
  # Save the cCRE ids corresponding to the consistent marking by the corresponding mark
  consistent[[i]] <- tmp[tmp$consecutive == 'consistent', ]$V1
  
}

# Subset the cases that have a mark consistently for the set of activating marks
m <- m[m$cCRE_id %in% Reduce(union, list(consistent$H3K27ac, 
                                         consistent$H3K9ac, 
                                         consistent$H3K4me3)), ]


# Add ATAC-seq data
c <- read.table('ATACseq/cCREs.ATACseq_overlap.tsv')

## Convert it to short format
c <- c[, c(6:8, 13)]
c$V13 <- 1
c <- reshape(unique(c), idvar = c('V6', 'V7'), timevar = "V8", direction = "wide")
c <- c[, -3]
c <- c[, c('V6', 'V7', 'V13.H000', 'V13.H012', 'V13.H024', 'V13.H096')]
c[is.na(c)] <- 0
rownames(c) <- c$V7
c <- c[, 3:6]
c[c > 0] <- 1 

## Obtain time point of first appearance
### Determine the equivalente between timepoints in ChIP-seq and ATAC-seq
tps <- c(1, 5, 7, 10.5)
c$appear <- apply(c, 1, function(a) tps[which(a!=0)[1]])
c[is.na(c)] <- 0
c$cCRE_id <- rownames(c)

m$supracategory <- 'CA/TF'
m[m$category == 'PLS', ]$supracategory <- 'PLS'
m[m$category == 'pELS', ]$supracategory <- 'pELS'
m[m$category == 'dELS', ]$supracategory <- 'dELS'

## Add class showing if the cCRE has ATAC-seq peak
c <- c[c$appear != 0, ]
m$ATAC <- 'closed'
m[m$cCRE_id %in% c$cCRE_id, ]$ATAC <- 'accessible'

# Set marks order and the geometric progression to encode each mark
ordered_marks <- c('H3K4me1', 'H3K4me2', 'H3K27ac', 'H3K9ac', 'H3K4me3')
code <- 2**(0:4)

tmp <- matrix(0, nrow = nrow(m), ncol = 12)
rownames(tmp) <- m$cCRE_id
tmp <- as.data.frame(tmp)

for ( i in 1:5 ) {
  n <- read.table(paste0(ordered_marks[i],
                         "/",
                         ordered_marks[i],
                         ".peaks.dynamics.binary.tsv"), stringsAsFactors = F)
  n[, 2:13] <- n[, 2:13]*code[i]
  rownames(n) <- n$V1
  n$V1 <- NULL
  n <- merge(tmp, n,
             by = 'row.names',
             all.x = TRUE)
  rownames(n) <- n$Row.names
  n <- n[, 14:25]
  n[is.na(n)] <- 0
  tmp <- tmp + n[rownames(tmp), ]
  tmp[is.na(tmp)] <- 0
}

tmp$cCRE_id <- rownames(tmp)
rownames(tmp) <- NULL
tmp <- tmp[, c(13, 1:12)]

# Check for minor combinations
filter <- tmp %>%
  group_by(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10, V11, V12) %>%
  dplyr::summarise(count = n()) %>%
  as.data.frame()

# Remove combinations that appear just once
filter <- merge(tmp, filter,
                all.x = TRUE)
# tmp <- tmp[tmp$cCRE_id %in% filter[filter$count > 1, ]$cCRE_id, ]
colnames(tmp)[2:13] <- c("H000", "H003", "H006",
                         "H009", "H012", "H018",
                         "H024", "H036", "H048",
                         "H072", "H120", "H168")

# Plot
## Obtain combinations and their associated identifier
class <- unlist(sapply(1:6, function(y) combn(ordered_marks, y, paste, collapse = "+")))
names(class) <- as.character(unlist(sapply(1:6, function(y) combn(c(1, 2, 4, 8, 16, 32), y, sum))))
class <- c(class, c(' 0' = 'unmarked'))

color <- rep('darkgrey', 31)

# Add custom colors based on the number of marks that are making up the specific combination
color <- ifelse(stringr::str_count(class, "\\+") == '0', yes = '#fcbba1',
                no = ifelse(stringr::str_count(class, "\\+") == '1', yes = '#fc9272',
                            no = ifelse(stringr::str_count(class, "\\+") == '2', yes = '#fb6a4a',
                                        no = ifelse(stringr::str_count(class, "\\+") == '3', yes = '#de2d26',
                                                    no = ifelse(stringr::str_count(class, "\\+") == '4', yes = '#a50f15',
                                                                no = ifelse(stringr::str_count(class, "\\+") == '5', yes = '#67000d',
                                                                            no = "#fee5d9"))))))

names(color) <- as.character(unlist(sapply(1:6, function(y) combn(c(1, 2, 4, 8, 16, 32), y, sum))))

color[64] <- "#fee5d9"
names(color)[64] <- '0'

names(color) <- str_pad(names(color), 2, 'left', ' ')
names(class) <- str_pad(names(class), 2, 'left', ' ')

# Generate label colors
lab_color <- color
lab_color <- gsub('#.*', 'black', lab_color)

# Add number identifier to the class labels
class <- setNames(paste0(class, ' (', str_replace_all(names(class), ' ', ''), ')'), 
                  nm = names(class))

# Generate information on marking by H3K9me3 and H3K27me3
repressive_marks <- c('H3K9me3', 'H3K27me3')

barplot_silent <- matrix()
barplot_silent <- matrix(0, nrow = nrow(m), ncol = 12)
rownames(barplot_silent) <- m$cCRE_id
barplot_silent <- as.data.frame(barplot_silent)

## Obtain combinations and their associated identifier
repressive_class <- unlist(sapply(1:2, function(y) combn(repressive_marks, y, paste, collapse = "+")))
names(repressive_class) <- as.character(unlist(sapply(1:2, function(y) combn(c(1, 2), y, sum))))
repressive_class <- c(repressive_class, c(' 0' = 'unmarked'))

repressive_color <- rep('darkgrey', 4)
names(repressive_color) <- as.character(unlist(sapply(1:2, function(y) combn(c(1, 2), y, sum))))

repressive_color['0'] <- "lightgrey"
repressive_color['1'] <- "#A7ADD4"
repressive_color['2'] <- "#1D2976"
repressive_color['3'] <- "black"

names(repressive_color) <- str_pad(names(repressive_color), 2, 'left', ' ')
names(repressive_class) <- str_pad(names(repressive_class), 2, 'left', ' ')

# Generate label colors
repressive_lab_color <- color
repressive_lab_color <- gsub('#.*', 'black', lab_color)

# Add number identifier to the class labels
repressive_class <- setNames(paste0(repressive_class, ' (', str_replace_all(names(repressive_class), ' ', ''), ')'), 
                             nm = names(repressive_class))

for ( i in 1:2 ) {
  n <- read.table(paste0(repressive_marks[i],
                         "/",
                         repressive_marks[i],
                         ".peaks.dynamics.binary.tsv"), stringsAsFactors = F)
  n[, 2:13] <- n[, 2:13]*code[i]
  rownames(n) <- n$V1
  n$V1 <- NULL
  n <- merge(m, n,
             by.x = 'cCRE_id',
             by.y = 'row.names',
             all.x = TRUE)
  rownames(n) <- n$cCRE_id
  n$cCRE_id <- NULL
  n <- n[, 13:24]
  n[is.na(n)] <- 0
  barplot_silent <- barplot_silent + n[rownames(barplot_silent), ]
  barplot_silent[is.na(barplot_silent)] <- 0
}

colnames(barplot_silent) <- c("H000", "H003", "H006",
                              "H009", "H012", "H018",
                              "H024", "H036", "H048",
                              "H072", "H120", "H168")

# Generate plot for accessible cCREs
plot <- do.call(rbind, apply(tmp[tmp$cCRE_id %in% m[m$ATAC == 'accessible', ]$cCRE_id, 
                                 c('cCRE_id', 'H000', 'H024', 'H048', 'H072', 'H120', 'H168')], 1, function(x) {
                                   x <- na.omit(x[-1])
                                   data.frame(x = names(x), node = x, 
                                              next_x = dplyr::lead(names(x)), 
                                              next_node = dplyr::lead(x), row.names = NULL)
                                 })) %>%
  mutate(x = factor(x, names(tmp[, c('cCRE_id', 'H000', 'H024', 'H048', 'H072', 'H120', 'H168')])[-1]),
         next_x = factor(next_x, names(tmp[, c('cCRE_id', 'H000', 'H024', 'H048', 'H072', 'H120', 'H168')])[-1]))

plot$node <- str_pad(plot$node, 2, 'left', ' ')
plot$next_node <- str_pad(plot$next_node, 2, 'left', ' ')

plot$label <- ''
plot$label <- repressive_class[plot$node]
plot$label <- gsub('\\+', '\n', plot$label)

plot$label <- gsub("\\s*\\([^\\)]+\\)","",plot$label)

# Add extended names for categories that represent up to the 95% of the total 
## Select this categories
freq <- as.data.frame(table(as.matrix(tmp[tmp$cCRE_id %in% m[m$ATAC == 'accessible', ]$cCRE_id, 
                                          2:13])))
freq$perc <- freq$Freq/sum(freq$Freq)*100
freq <- freq[order(freq$perc, decreasing = TRUE), ]
freq$cumsum <- cumsum(freq$perc)

selected_combs <- freq[1:(nrow(freq[freq$cumsum < 95, ]) + 1), ]$Var1
selected_combs <- str_pad(selected_combs, 2, 'left', ' ')

# Label only the most common combinations
plot$label <- ''
plot[plot$node %in% selected_combs, ]$label <- class[plot[plot$node %in% selected_combs, ]$node]
plot[plot$node %in% selected_combs, ]$label <- gsub('\\+', '\n', plot[plot$node %in% selected_combs, ]$label)

plot$label <- gsub("\\s*\\([^\\)]+\\)","",plot$label)

new_lab_color <- lab_color
new_lab_color[] <- 'darkgrey'
new_lab_color[names(new_lab_color) %in% selected_combs] <- 'black'

pdf('fig_S2cd.pdf',
    height = 4,
    width = 6)

print(ggplot(data = plot,
             aes(x = x,
                 next_x = next_x,
                 node = node,
                 next_node = next_node,
                 label = label,
                 fill = node)) +
        geom_sankey(flow.alpha = 0.5,
                    # node.color = NA,
                    node.color = 'black',
                    show.legend = TRUE) +
        ggtitle(paste0(format(nrow(tmp[tmp$cCRE_id %in% m[m$ATAC == 'accessible', ]$cCRE_id, ]),
                              big.mark = ",",
                              scientific=FALSE),
                       ' accessible cCREs')) +
        geom_sankey_text(aes(color = node),
                         size = 1.5,
                         hjust = 0,
                         lineheight = 0.7,
                         position = position_nudge(x = 0.1)) +
        scale_color_manual(values = new_lab_color) +
        scale_fill_manual(labels = class,
                          values = color) +
        guides(fill = guide_legend(ncol = 1),
               color = 'none') +
        xlab('') +
        theme_sankey() +
        theme(axis.text.x = element_text(angle = 90, 
                                         vjust = 0.5, 
                                         hjust = 1),
              legend.position = 'none',
              legend.text = element_text(size = 5),
              legend.key.size = unit(3, "mm"),
              plot.title = element_text(hjust = 0.5),
              plot.margin = margin(2, 30, 2, 30, "mm"))
)

## Add data on overlap from H3K9me3 and H3K27me3
tmp.barplot_silent <- barplot_silent[m[m$ATAC == 'accessible', ]$cCRE_id, 
                                     c('H000', 'H024', 'H048', 'H072', 'H120', 'H168')] %>%
  melt() %>%
  group_by_all() %>%
  dplyr::summarise(count = n()) %>%
  as.data.frame()
tmp.barplot_silent$value <- factor(tmp.barplot_silent$value,
                                   levels = c('0', '1', '2', '3'))

# Do alluvial plot of repressive marks
plot <- do.call(rbind, apply(barplot_silent[m[m$ATAC == 'accessible', ]$cCRE_id, 
                                            c('H000', 'H024', 'H048', 'H072', 'H120', 'H168')], 1, function(x) {
                                              data.frame(x = names(x), node = x, 
                                                         next_x = dplyr::lead(names(x)), 
                                                         next_node = dplyr::lead(x), row.names = NULL)
                                            })) %>%
  mutate(x = factor(x, names(tmp[, c('H000', 'H024', 'H048', 'H072', 'H120', 'H168')])),
         next_x = factor(next_x, names(tmp[, c('H000', 'H024', 'H048', 'H072', 'H120', 'H168')])))

plot$node <- str_pad(plot$node, 2, 'left', ' ')
plot$next_node <- str_pad(plot$next_node, 2, 'left', ' ')

plot$label <- ''
plot$label <- repressive_class[plot$node]
plot$label <- gsub('\\+', '\n', plot$label)

plot$label <- gsub("\\s*\\([^\\)]+\\)","",plot$label)

repressive_new_lab_color <- repressive_lab_color
repressive_new_lab_color[] <- 'black'

# Generate the plot only for cases marked by silent marks in at least one time point
barplot_silent_ss <- barplot_silent[apply(barplot_silent, 1, function(a) sum(a) > 0), ]
plot <- do.call(rbind, apply(barplot_silent_ss[intersect(m[m$ATAC == 'accessible', ]$cCRE_id, 
                                                         rownames(barplot_silent_ss)), 
                                               c('H000', 'H024', 'H048', 'H072', 'H120', 'H168')], 1, function(x) {
                                                 data.frame(x = names(x), node = x, 
                                                            next_x = dplyr::lead(names(x)), 
                                                            next_node = dplyr::lead(x), row.names = NULL)
                                               })) %>%
  mutate(x = factor(x, names(tmp[, c('H000', 'H024', 'H048', 'H072', 'H120', 'H168')])),
         next_x = factor(next_x, names(tmp[, c('H000', 'H024', 'H048', 'H072', 'H120', 'H168')])))

plot$node <- str_pad(plot$node, 2, 'left', ' ')
plot$next_node <- str_pad(plot$next_node, 2, 'left', ' ')

plot$label <- ''
plot$label <- repressive_class[plot$node]
plot$label <- gsub('\\+', '\n', plot$label)

plot$label <- gsub("\\s*\\([^\\)]+\\)","",plot$label)

repressive_new_lab_color <- repressive_lab_color
repressive_new_lab_color[] <- 'black'

print(ggplot(data = plot,
             aes(x = x,
                 next_x = next_x,
                 node = node,
                 next_node = next_node,
                 label = label,
                 fill = node)) +
        geom_sankey(flow.alpha = 0.5,
                    # node.color = NA,
                    node.color = 'black',
                    show.legend = TRUE) +
        ggtitle(paste0(format(nrow(barplot_silent_ss[intersect(m[m$ATAC == 'accessible', ]$cCRE_id, 
                                                               rownames(barplot_silent_ss)), ]),
                              big.mark = ",",
                              scientific=FALSE),
                       ' accessible cCREs')) +
        geom_sankey_text(aes(color = node),
                         size = 1.5,
                         hjust = 0,
                         lineheight = 0.7,
                         position = position_nudge(x = 0.1)) +
        scale_color_manual(values = repressive_new_lab_color) +
        scale_fill_manual(labels = repressive_class,
                          values = repressive_color) +
        guides(fill = guide_legend(ncol = 1),
               color = 'none') +
        xlab('') +
        theme_sankey() +
        theme(axis.text.x = element_text(angle = 90, 
                                         vjust = 0.5, 
                                         hjust = 1),
              legend.position = 'none',
              legend.text = element_text(size = 5),
              legend.key.size = unit(3, "mm"),
              plot.title = element_text(hjust = 0.5),
              plot.margin = margin(2, 30, 2, 30, "mm"))
)
dev.off()
