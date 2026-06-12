#!/usr/bin/env r

#******************************************************
## Test the causal effect that gene activation has on
## the signal of different histone PTMs at promoters
#******************************************************


#-------------------------
## Load required libraries
#-------------------------
library(tidyverse)
library(vroom)
library(tjbal)


#------------------
## Helper functions
#------------------

# Check if a variable is monotonically increasing (or decreasing)
mon <-
  function(x, type) {
    if (type == "inc") {
      return(all(x == cummax(x)))
    }
    if (type == "dec") {
      return(all(x == cummin(x)))
    }
  }

#-----------
## Load Data
#-----------

# RNA Expression Table
gene_exp <-
  vroom("data/selected.genes.rep.2.3.after.QN.merged.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    )
  )

# List of genes that get activated during transdifferentiation
activated_gene_set <-
  vroom("data/257.notExpressed.0h.txt",
    col_names = "gene_id",
    delim = "\t"
  )

# ChIP-seq signal for different histone PTMs
h4k20me1 <-
  vroom("data/updated_peaks/H4K20me1.matrix.after.QN.merged.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    )
  )

h3k4me1 <-
  vroom("data/updated_peaks/H3K4me1.matrix.after.QN.merged.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    )
  )

h3k4me2 <-
  vroom("data/updated_peaks/H3K4me2.matrix.after.QN.merged.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    )
  )

h3k4me3 <-
  vroom("data/updated_peaks/H3K4me3.matrix.after.QN.merged.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    )
  )

h3k27ac <-
  vroom("data/updated_peaks/H3K27ac.matrix.after.QN.merged.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    )
  )

h3k36me3 <-
  vroom("data/updated_peaks/H3K36me3.matrix.after.QN.merged.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    )
  )

h3k27me3 <-
  vroom("data/updated_peaks/H3K27me3.matrix.after.QN.merged.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    )
  )

h3k9me3 <-
  vroom("data/updated_peaks/H3K9me3.matrix.after.QN.merged.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    )
  )

h3k9ac <-
  vroom("data/updated_peaks/H3K9ac.matrix.after.QN.merged.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    )
  )


#------------------
## Pre-process Data
#------------------

# For every histone PTM filter for genes present
# in the gene_exp table and change the fromat of
# the table from "wider" to "longer"
h3k27ac <-
  h3k27ac %>%
  pivot_longer(-gene_id, names_to = "time", values_to = "H3K27ac") %>%
  filter(gene_id %in% gene_exp$gene_id) %>%
  mutate(time = as.integer(time))

h3k4me1 <-
  h3k4me1 %>%
  pivot_longer(-gene_id, names_to = "time", values_to = "H3K4me1") %>%
  filter(gene_id %in% gene_exp$gene_id) %>%
  mutate(time = as.integer(time))

h3k4me2 <-
  h3k4me2 %>%
  pivot_longer(-gene_id, names_to = "time", values_to = "H3K4me2") %>%
  filter(gene_id %in% gene_exp$gene_id) %>%
  mutate(time = as.integer(time))

h3k4me3 <-
  h3k4me3 %>%
  pivot_longer(-gene_id, names_to = "time", values_to = "H3K4me3") %>%
  filter(gene_id %in% gene_exp$gene_id) %>%
  mutate(time = as.integer(time))

h3k36me3 <-
  h3k36me3 %>%
  pivot_longer(-gene_id, names_to = "time", values_to = "H3K36me3") %>%
  filter(gene_id %in% gene_exp$gene_id) %>%
  mutate(time = as.integer(time))

h4k20me1 <-
  h4k20me1 %>%
  pivot_longer(-gene_id, names_to = "time", values_to = "H4K20me1") %>%
  filter(gene_id %in% gene_exp$gene_id) %>%
  mutate(time = as.integer(time))

h3k9ac <-
  h3k9ac %>%
  pivot_longer(-gene_id, names_to = "time", values_to = "H3K9ac") %>%
  filter(gene_id %in% gene_exp$gene_id) %>%
  mutate(time = as.integer(time))

h3k9me3 <-
  h3k9me3 %>%
  pivot_longer(-gene_id, names_to = "time", values_to = "H3K9me3") %>%
  filter(gene_id %in% gene_exp$gene_id) %>%
  mutate(time = as.integer(time))

h3k27me3 <-
  h3k27me3 %>%
  pivot_longer(-gene_id, names_to = "time", values_to = "H3K27me3") %>%
  filter(gene_id %in% gene_exp$gene_id) %>%
  mutate(time = as.integer(time))

# Merge data for all histone PTMs
marks <-
  h3k4me1 %>%
  full_join(h3k4me2) %>%
  full_join(h3k4me3) %>%
  full_join(h3k27ac) %>%
  full_join(h3k36me3) %>%
  full_join(h4k20me1) %>%
  full_join(h3k9ac) %>%
  full_join(h3k27me3) %>%
  full_join(h3k9me3) %>%
  pivot_longer(
    -c("gene_id", "time"),
    names_to = "mark", values_to = "signal"
  ) %>%
  mutate(mark = factor(mark, levels = mark %>% unique))

# Check the distribution of signal per mark
marks %>%
  mutate(signal = log2(1 + signal)) %>%
  group_by(mark) %>%
  summarise(
    min = quantile(signal, probs = 0),
    q25 = quantile(signal, probs = 0.25),
    med = quantile(signal, probs = 0.5),
    q75 = quantile(signal, probs = 0.75),
    max = quantile(signal, probs = 1),
  )

# Set a threshold for RNA expression to define gene activation
thresh <- 1

# Reshape gene expression table
gene_exp <-
  gene_exp %>%
  pivot_longer(-gene_id, names_to = "time", values_to = "gene_exp") %>%
  mutate(
    time = as.integer(time),
    gene_exp_thresh = (gene_exp > thresh) %>% as.integer
  )

# Distribution of the number of timepoints that genes pass the cutoff
gene_exp %>%
  group_by(gene_id) %>%
  summarise(sum = sum(gene_exp_thresh)) %>%
  pull(sum) %>%
  table()


#------------------------------
## Set control and tested genes
#------------------------------

# Find genes that are active in all timepoints
control4disappearance <-
  gene_exp %>%
  group_by(gene_id) %>%
  summarise(sum = sum(gene_exp_thresh, na.rm = TRUE)) %>%
  filter(sum == 12) %>%
  select(-sum) %>%
  ungroup()

# Find genes that are not active in any timepoint
# These are the silent genes
control4appearance <-
  gene_exp %>%
  group_by(gene_id) %>%
  summarise(sum = sum(gene_exp_thresh, na.rm = TRUE)) %>%
  filter(sum == 0) %>%
  select(-sum) %>%
  ungroup()


# Find out genes that get activated at some timepoint
# and remain active after that
test4appearance <-
  gene_exp %>%
  filter(!(
    gene_id %in% c(control4appearance$gene_id, control4disappearance$gene_id)
  )) %>%
  group_by(gene_id) %>%
  summarise(keep = mon(gene_exp_thresh, "inc")) %>%
  filter(keep) %>%
  select(-keep) %>%
  ungroup()

# Out of the set of 257 activated genes, 20 are not included in the tested
# genes. That is because they pass the cutoff 1 TPM at a certain timepoint, but
# then they go bellow it at a later timepoint.
#
# Fix this by considering only the first time they pass the cutoff.
# Then, redefine the tested genes to include only this set of activated genes.
gene_exp <-
  gene_exp %>%
  mutate(
    to_fix = gene_id %in%
      setdiff(activated_gene_set$gene_id, test4appearance$gene_id)
  ) %>%
  group_by(gene_id) %>%
  mutate(cumany = cumany(gene_exp_thresh)) %>%
  ungroup() %>%
  mutate(
    gene_exp_thresh = ifelse(to_fix & cumany & gene_exp_thresh == 0, 1, gene_exp_thresh)
  ) %>%
  select(gene_id:gene_exp_thresh)

test4appearance <-
  gene_exp %>%
  filter(!(
    gene_id %in% c(control4appearance$gene_id, control4disappearance$gene_id)
  )) %>%
  group_by(gene_id) %>%
  summarise(keep = mon(gene_exp_thresh, "inc")) %>%
  filter(keep & gene_id %in% activated_gene_set$gene_id) %>%
  select(-keep)


# Join expression and histone marks table
data <-
  gene_exp %>%
  right_join(marks)

# Data to test
data_appearance <-
  data %>%
  right_join(
    rbind(control4appearance, test4appearance)
  )


#--------------------
## Causality Analysis
#--------------------
out_causality_appearance <- list()
trim_npre <- 2

for (m in unique(data_appearance$mark)) {
  print(m)

  out_causality_appearance[[m]] <-
    tjbal(
      signal ~ gene_exp_thresh,
      data = data_appearance %>% filter(mark == m),
      index = c("gene_id", "time"), estimator = "kernel",
      demean = FALSE, vce = "jackknife", nsims = 1000,
      cores = 1, parallel = FALSE, seed = 1234, trim.npre = trim_npre
    )

}
