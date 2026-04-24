# Only on marks preceding gene expression ----
set.seed(123)
setwd("../analysis/Random-Forest/H3K4me12_H3K27ac/")

# Performance data (upreg genes, test set)
perf <- fread("summary.txt", data.table = F, fill = TRUE)
perf <- perf[, 1:4]
tps <- c("H000","H003","H006","H009","H012","H018","H024","H036","H048") 

# Subset timepoints, set as numeric
perf <- subset(perf, tc %in% tps & te %in% tps)  
perf$te <- as.numeric(gsub("H", "", perf$te))         
perf$tc <- as.numeric(gsub("H", "", perf$tc))

# Define groups: backward, forward, same
perf$timing <- NA
for (i in 1:nrow(perf)){
  perf[i, "timing"] = ifelse(perf[i, "tc"] > perf[i, "te"], "Backward", "Forward")
  if(perf[i, "tc"] == perf[i, "te"]){perf[i, "timing"] = "Same"}
}

# Mean RMSE per group
perf2 <- perf %>% 
  group_by(te, timing) %>% 
  summarise(mean = mean(rmse)) %>%
  as.data.frame()

# Compute ratios and build data.frame to plot
Forward <- c(subset(perf2, timing == "Forward")$mean, NA)/subset(perf2, timing == "Same")$mean   # NA added (48h has not forward)
Backward <- c(NA, subset(perf2, timing == "Backward")$mean)/subset(perf2, timing == "Same")$mean  # NA added (0h has not backward)
Same <- rep(1, 9)

perf3 <- data.frame(Backward, Forward, te = unique(perf2$te))
perf3m <- melt(perf3, id.vars = "te")
perf3m <- perf3m[!is.na(perf3m$value), ] # Remove 0h backward and 48h forward


# Plot
p  <- ggplot(perf3m, aes(variable, value, label = as.character(te))) +
  geom_boxplot(fill = "gray95") +
  #geom_hline(yintercept = 1, col = 2, lty = 2) +
  geom_signif(comparisons = list(c("Forward", "Backward")), 
              tip_length = 0.01, y_position = 1.07, map_signif_level = F, 
              annotations = sprintf("p = %0.2e", wilcox.test(perf3m$value ~ perf3m$variable)$p.value), textsize = 6) +
  geom_text_repel(seed = 123, size = 6) +
  theme_bw(base_size = 24) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_blank()
  ) +
  ylim(0.995, 1.075) +
  labs(y = expression(bar("RMSE")*"/RMSE"["(i=j)"])) 

pdf("fig_3f.pdf", 
    width = 6.25, 
    height = 4.75)
p

# Only on marks delayed to gene expression ----
set.seed(123)
setwd("../analysis/Random-Forest/H3K9ac_H3K4me3_H3K36me3_H4K20me1/")

# Performance data (upreg genes, test set)
perf <- fread("summary.txt", data.table = F, fill = TRUE)
perf <- perf[, 1:4]
tps <- c("H000","H003","H006","H009","H012","H018","H024","H036","H048") 

# Subset timepoints, set as numeric
perf <- subset(perf, tc %in% tps & te %in% tps)  
perf$te <- as.numeric(gsub("H", "", perf$te))         
perf$tc <- as.numeric(gsub("H", "", perf$tc))

# Define groups: backward, forward, same
perf$timing <- NA
for (i in 1:nrow(perf)){
  perf[i, "timing"] = ifelse(perf[i, "tc"] > perf[i, "te"], "Backward", "Forward")
  if(perf[i, "tc"] == perf[i, "te"]){perf[i, "timing"] = "Same"}
}

# Mean RMSE per group
perf2 <- perf %>% 
  group_by(te, timing) %>% 
  summarise(mean = mean(rmse)) %>%
  as.data.frame()

# Compute ratios and build data.frame to plot
Forward <- c(subset(perf2, timing == "Forward")$mean, NA)/subset(perf2, timing == "Same")$mean   # NA added (48h has not forward)
Backward <- c(NA, subset(perf2, timing == "Backward")$mean)/subset(perf2, timing == "Same")$mean  # NA added (0h has not backward)
Same <- rep(1, 9)

perf3 <- data.frame(Backward, Forward, te = unique(perf2$te))
perf3m <- melt(perf3, id.vars = "te")
perf3m <- perf3m[!is.na(perf3m$value), ] # Remove 0h backward and 48h forward


# Plot
p  <- ggplot(perf3m, aes(variable, value, label = as.character(te))) +
  geom_boxplot(fill = "gray95") +
  #geom_hline(yintercept = 1, col = 2, lty = 2) +
  geom_signif(comparisons = list(c("Forward", "Backward")), 
              tip_length = 0.01, y_position = 1.14, map_signif_level = F, 
              annotations = sprintf("p = %0.2e", wilcox.test(perf3m$value ~ perf3m$variable)$p.value), textsize = 6) +
  geom_text_repel(seed = 123, size = 6) +
  theme_bw(base_size = 24) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"),
        axis.text = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.title.x = element_blank()
  ) +
  ylim(0.935, 1.16) +
  labs(y = expression(bar("RMSE")*"/RMSE"["(i=j)"])) 

p

dev.off()
