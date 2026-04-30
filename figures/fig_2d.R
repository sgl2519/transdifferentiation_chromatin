library(ggplot2)
library(dplyr)

# Load matrix of pairs and their profile
df <- read.table('../analysis/cCREs/accessible.pairs_wrt.HiC.support_get20kb.let100kb_profiles.tsv.gz',
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

# Compute profile sharing on cCRE pairs with constant support or unsupport ----
df1 <- df
df1$correspondence <- 'different'
df1[df1$pf_A == df1$pf_B, ]$correspondence <- 'same'

df1 <- df1 %>%
  group_by(cCRE_B, cCRE_A, pair, dist.clusters, correspondence) %>%
  dplyr::summarise(count = n()) %>%
  as.data.frame()

# Extract stably-contacting cases
df1 <- df1[df1$count == 8, ]

# Summarise
df1 <- df1 %>%
  group_by(pair, dist.clusters, correspondence) %>%
  dplyr::summarise(count = n()) %>%
  as.data.frame()

# Compute proportions of shared profile
df1 <- df1 %>%
  group_by(pair, dist.clusters) %>%
  dplyr::summarise(prop = count/sum(count),
                   count = count,
                   correspondence = correspondence) %>%
  as.data.frame()

df1$pair <- factor(df1$pair,
                   levels = c('no_HiC.support', 'HiC.support'))


# Compute paired t-test
test <- t.test(df1[df1$correspondence == 'same' &
                     df1$pair == 'HiC.support', ]$prop, 
               df1[df1$correspondence == 'same' &
                     df1$pair == 'no_HiC.support', ]$prop, 
               # paired = TRUE,
               alternative = 'greater')$p.value

test <- formatC(test, format = "e", digits = 2)
pdf('fig_2d.pdf',
    height = 6,
    width = 9.5)
ggplot(data = df1[df1$correspondence == 'same', ],
       aes(x = prop,
           y = pair,
           color = pair)) +
  geom_violin(fill = 'white') +
  # scale_fill_manual(values = c('HiC.support' = '#2ca25f',
  #                              'no_HiC.support' = '#bdbdbd')) +
  geom_boxplot(fill = 'white',
               width = 0.3) +
  scale_color_manual(values = c('HiC.support' = '#2ca25f',
                                'no_HiC.support' = '#bdbdbd'),
                     labels = c('no_HiC.support' = 'Hi-C-unsupported\npairs',
                                'HiC.support' = 'Hi-C-supported\npairs')) +
  geom_jitter(data = df1[df1$correspondence == 'same', ],
             aes(x = prop,
                 y = pair,
                 fill = dist.clusters),
             shape = 21,
             color = 'black',
             height = 0.1) +
  geom_text(data = as.data.frame(test),
            label = paste0('p-value = ', test),
            color = 'black',
            x = 0.15,
            y = 1.5,
            hjust = 0,
            size = 3) +
  scale_fill_manual(values = c("[20, 30 kb)" = '#023858',
                               "[30, 40 kb)" = '#045a8d',
                               "[40, 50 kb)" = '#0570b0', 
                               "[50, 60 kb)" = '#3690c0',
                               "[60, 70 kb)" = '#74a9cf',
                               "[70, 80 kb)" = '#a6bddb',
                               "[80, 90 kb)" = '#d0d1e6',
                               "[90, 100 kb]" = '#ece7f2')) +
  scale_x_continuous(name = 'Proportion of cCRE pairs\nsharing chromatin trajectory',
                     labels = scales::percent) +
  scale_y_discrete(name = '',
                   labels = c('no_HiC.support' = 'Hi-C-unsupported\npairs',
                              'HiC.support' = 'Hi-C-supported\npairs')) +
  guides(color = guide_legend(title = "Chromatin trajectory\ncorrespondence\nwithin cCRE pairs"),
         fill = guide_legend(title = "Distance between\ncCREs")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, 
                                   hjust = 1),
        strip.background = element_rect(fill = "white"),
        plot.margin = unit(c(2, 2, 2, 2), "inch"))
