# Represent the proportion of shared profiles based on the distance cCREs are separated by, disregarding Hi-C contact ----

# Read table of cCREs pairs
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

df1 <- df
df1$correspondence <- 'different'
df1[df1$pf_A == df1$pf_B, ]$correspondence <- 'same'
df1$pf <- 'different'
df1[df1$correspondence == 'same', ]$pf <- df1[df1$correspondence == 'same', ]$pf_A

# Compute number of cCRE pairs based on the profile they share/not share
count.label <- c(table(df1$pf)/8)

# Plot
pdf('fig_S5a.pdf',
    height = 6,
    width = 10)

# Compute number of cases that share or not the same profile in the vicinity of cCREs
plot <- df1 %>%
  group_by(dist.clusters, correspondence) %>%
  dplyr::summarise(count = n()/8) %>%
  group_by(dist.clusters) %>%
  dplyr::summarise(prop = count/sum(count),
                   total = sum(count),
                   correspondence = correspondence) %>%
  as.data.frame()

# Plot
ggplot(data = plot[plot$correspondence == 'same', ],
       aes(y = dist.clusters,
           x = prop,
           fill = dist.clusters)) +
  geom_col(color = 'black',
           alpha = 0.75) +
  scale_x_continuous(name = 'Proportion of cCRE pairs sharing chromatin trajectory',
                     labels = scales::percent) +
  ylab('Distance between cCREs') +
  scale_fill_manual(name = 'Distance between\ncCREs',
                    values = c("[20, 30 kb)" = '#023858',
                               "[30, 40 kb)" = '#045a8d',
                               "[40, 50 kb)" = '#0570b0', 
                               "[50, 60 kb)" = '#3690c0',
                               "[60, 70 kb)" = '#74a9cf',
                               "[70, 80 kb)" = '#a6bddb',
                               "[80, 90 kb)" = '#d0d1e6',
                               "[90, 100 kb]" = '#ece7f2')) +
  geom_text(aes(y = dist.clusters,
                label = paste0(gsub(' ',
                                    '',
                                    format(total, 
                                           big.mark = ",")), 
                               ' cCRE pairs')),
            x = max(plot[plot$correspondence == 'same', ]$prop) + 0.02, # Set text's position to the right end of the plot
            color = 'black',
            hjust = 0,
            size = 3,
            fontface = 'italic') +
  coord_cartesian(clip = 'off') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 1),
        legend.position = 'none',
        plot.margin = unit(c(2, 3.5, 2, 3.5), "inch"))

dev.off()
