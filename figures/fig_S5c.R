library(ggplot2)
library(dplyr)

# Load matrix of pairs and their profile
df <- read.table('../analysis/cCREs/tHMM/accessible.pairs_wrt.HiC.support_get20kb.let100kb_profiles.tsv.gz',
                 header = TRUE)

# Add column indicating distance between cCREs
df$distance <- df$start_B - df$end_A

# Stratify distance by sections
df$dist.clusters <- 'none'
df[(df$distance >= 20000) &
     (df$distance < 30000), ]$dist.clusters <- '[20, 30 kb)'
df[(df$distance >= 30000) &
     (df$distance < 40000), ]$dist.clusters <- '[30, 40 kb)'
df[(df$distance >= 40000) &
     (df$distance < 50000), ]$dist.clusters <- '[40, 50 kb)'
df[(df$distance >=50000) &
     (df$distance < 60000), ]$dist.clusters <- '[50, 60 kb)'
df[(df$distance >= 60000) &
     (df$distance < 70000), ]$dist.clusters <- '[60, 70 kb)'
df[(df$distance >= 70000) &
     (df$distance < 80000), ]$dist.clusters <- '[70, 80 kb)'
df[(df$distance >= 80000) &
     (df$distance < 90000), ]$dist.clusters <- '[80, 90 kb)'
df[(df$distance >= 90000) &
     (df$distance <= 100000), ]$dist.clusters <- '[90, 100 kb]'

# Obtain number of pairs separated by each distance
label <- as.data.frame(table(df$dist.clusters)/8)
label <- setNames(paste0(label$Var1, '\n', 
                         gsub(' ', '', 
                              format(label$Freq, 
                                     big.mark   = ',')), 
                         ' pairs'),
                  nm = label$Var1)


# Plot distance distribution
pdf('fig_S5c.pdf',
    height = 6,
    width = 9.5)

# Compute proportion of shared profiles as a function of time points in contact ----
df1 <- df
df1$correspondence <- 'different'
df1[df1$pf_A == df1$pf_B, ]$correspondence <- 'same'

# Classify cCRE pairs based on the number of time points during which they are in contact
df1 <- df1 %>%
  group_by(cCRE_B, cCRE_A, pair, dist.clusters, correspondence) %>%
  dplyr::summarise(count = n()) %>%
  as.data.frame()

df1 <- df1[(df1$pair == 'HiC.support' & df1$count %in% 1:8) |
             (df1$pair == 'no_HiC.support' & df1$count == 8), ]

# Generate new column with grouping based on number of time points in contact
df1$group <- paste0(df1$pair, '_', df1$count)
df1 <- df1 %>%
  group_by(dist.clusters, correspondence, group) %>%
  dplyr::summarise(count = n()) %>%
  as.data.frame()

# Compute proportion of shared profiles per group and distance
plot <- df1 %>%
  group_by(dist.clusters, group) %>%
  dplyr::summarise(prop = count/sum(count),
                   total = sum(count),
                   correspondence = correspondence) %>%
  as.data.frame()

# Generate table with total number of contacts in each category
label <- plot %>%
  group_by(group) %>%
  dplyr::summarise(num = sum(total)/2) %>%
  as.data.frame()

# Add data for all cCRE pairs
all <- df1 %>%
  group_by(dist.clusters, correspondence) %>%
  dplyr::summarise(group = 'all',
                   count = sum(count)) %>%
  group_by(dist.clusters) %>%
  dplyr::summarise(group = 'all',
                   prop = count/sum(count),
                   total = sum(count),
                   correspondence = correspondence) %>%
  unique() %>%
  as.data.frame()

plot <- rbind(plot,
              all)

# Parse for plotting
plot$group <- factor(plot$group,
                     levels = c("all",
                                "no_HiC.support_8",
                                "HiC.support_1",
                                "HiC.support_2",
                                "HiC.support_3",
                                "HiC.support_4", 
                                "HiC.support_5",
                                "HiC.support_6",
                                "HiC.support_7",
                                "HiC.support_8"))
plot <- plot[order(plot$group), ]

## Plot
ggplot(data = plot[plot$correspondence == 'same', ],
       aes(x = prop,
           y = group,
           color = group,
           fill = group)) +
  geom_violin(color = 'black',
              alpha = 0.7) +
  geom_boxplot(width = 0.1,
               fill = 'white',
               color = '#525252') +
  xlab('Proportion of cCRE pairs sharing profile') +
  scale_y_discrete(labels = c("all" = "all cCRE pairs",
                              "no_HiC.support_8" = "contact in 0\ntime points",
                              "HiC.support_1" = "contact in 1\ntime point",
                              "HiC.support_2" = "contact in 2\ntime points",
                              "HiC.support_3" = "contact in 3\ntime points",
                              "HiC.support_4" = "contact in 4\ntime points", 
                              "HiC.support_5" = "contact in 5\ntime points",
                              "HiC.support_6" = "contact in 6\ntime points",
                              "HiC.support_7" = "contact in 7\ntime points",
                              "HiC.support_8" = "contact in 8\ntime points")) +
  scale_color_manual(values = c("all" = "grey",
                                "no_HiC.support_8" = "#ffffff",
                                "HiC.support_1" = "#ffffe5",
                                "HiC.support_2" = "#fff7bc",
                                "HiC.support_3" = "#fee391",
                                "HiC.support_4" = "#fec44f", 
                                "HiC.support_5" = "#fe9929",
                                "HiC.support_6" = "#ec7014",
                                "HiC.support_7" = "#cc4c02",
                                "HiC.support_8" = "#993404")) +
  scale_fill_manual(values = c("all" = "grey",
                               "no_HiC.support_8" = "#ffffff",
                               "HiC.support_1" = "#ffffe5",
                               "HiC.support_2" = "#fff7bc",
                               "HiC.support_3" = "#fee391",
                               "HiC.support_4" = "#fec44f", 
                               "HiC.support_5" = "#fe9929",
                               "HiC.support_6" = "#ec7014",
                               "HiC.support_7" = "#cc4c02",
                               "HiC.support_8" = "#993404")) +
  geom_jitter(data = plot[plot$correspondence == 'same', ],
             aes(x = prop,
                 y = group,
                 color = dist.clusters),
             height = 0.1, 
             size = 1) +
  scale_color_manual(values = c("[20, 30 kb)" = '#023858',
                                "[30, 40 kb)" = '#045a8d',
                                "[40, 50 kb)" = '#0570b0', 
                                "[50, 60 kb)" = '#3690c0',
                                "[60, 70 kb)" = '#74a9cf',
                                "[70, 80 kb)" = '#a6bddb',
                                "[80, 90 kb)" = '#d0d1e6',
                                "[90, 100 kb]" = '#ece7f2')) +
  ylab('') +
  theme_bw() +
  theme(legend.position = 'none',
        plot.margin = unit(c(1.5, 3, 1.5, 3), "inch"))

  
# Generate plot with proportion of chromatin profiles in the shared cases
df1 <- df
df1$correspondence <- 'different'
df1[df1$pf_A == df1$pf_B, ]$correspondence <- 'same'

# Classify cCRE pairs based on the number of time points during which they are in contact
df1 <- df1 %>%
  group_by(cCRE_B, cCRE_A, pair, dist.clusters, pf_A, pf_B) %>%
  dplyr::summarise(count = n()) %>%
  as.data.frame()

df1 <- df1[(df1$pair == 'HiC.support' & df1$count %in% 1:8) |
             (df1$pair == 'no_HiC.support' & df1$count == 8), ]

# Generate new column with grouping based on number of time points in contact
df1$group <- paste0(df1$pair, '_', df1$count)
df1$shared <- 'different'
df1[df1$pf_A == df1$pf_B, ]$shared <- df1[df1$pf_A == df1$pf_B, ]$pf_A

df1 <- df1 %>%
  group_by(shared, group) %>%
  dplyr::summarise(count = n()) %>%
  as.data.frame()

df1$group <- factor(df1$group,
                    levels = c("all",
                               "no_HiC.support_8",
                               "HiC.support_1",
                               "HiC.support_2",
                               "HiC.support_3",
                               "HiC.support_4", 
                               "HiC.support_5",
                               "HiC.support_6",
                               "HiC.support_7",
                               "HiC.support_8"))

df1$shared <- factor(df1$shared,
                     levels = rev(c('different', 8, 10, 9, 11, 12, 13,
                                    1, 2, 5, 6, 7, 3, 4)))

ggplot(data = df1[df1$shared != 'different', ],
       aes(x = count,
           y = group,
           fill = shared)) +
  geom_bar(position = "fill", 
           stat = "identity",
           alpha = 0.75,
           color = 'black') +
  scale_fill_manual(values = cols <- c("8" = "#bdbdbdff", 
                                       "10" = "black", 
                                       "9" = "#f7d68fff", 
                                       "11" = "#E6AB02", 
                                       "12" = "#e67f30ff", 
                                       "13" = "#E7298A",
                                       "1" = "#b2182b", 
                                       "2" = "#ef8a62", 
                                       "5" = "#fddbc7", 
                                       "6" = "#762a83",
                                       "7" = "#d1e5f0", 
                                       "3" = "#67a9cf", 
                                       "4" = "#2166ac")) +
  scale_y_discrete(name = '',
                   labels = c("all" = "all cCRE pairs",
                              "no_HiC.support_8" = "no contact",
                              "HiC.support_1" = "contact in 1\ntime point",
                              "HiC.support_2" = "contact in 2\ntime points",
                              "HiC.support_3" = "contact in 3\ntime points",
                              "HiC.support_4" = "contact in 4\ntime points", 
                              "HiC.support_5" = "contact in 5\ntime points",
                              "HiC.support_6" = "contact in 6\ntime points",
                              "HiC.support_7" = "contact in 7\ntime points",
                              "HiC.support_8" = "contact in 8\ntime points")) +
  scale_x_continuous(name = 'Proportion of shared\nchromatin trajectories',
                     labels = scales::percent) +
  theme_bw() +
  theme(legend.position = 'none',
        plot.margin = unit(c(1.5, 3, 1.5, 3), "inch"))

dev.off()
