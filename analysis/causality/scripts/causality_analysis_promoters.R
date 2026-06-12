#!/usr/bin/env r

#******************************************************
## Test the causal effect that a chromatin state change
## at promoters has on RNA expression
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
# It's used to identify the cases where there is a single chromatin state
# transition
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

# RNA Expression Table that includes lowly expressed genes
gene_exp_requantile <-
  vroom(
    "data/expression.matrix.quantile.renormalized.tsv",
    col_types = "cid"
  )


# List of silent genes that are not expressed
silenced_genes <-
  gene_exp %>%
  tail(n = 1552) %>%
  pull(gene_id)

# List of genes that are differentially expressed / changing
# during transdifferentiation
de_genes <-
  vroom("data/expression.matrix.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "000", "003", "006", "009", "012", "018",
      "024", "036", "048", "072", "120", "168"
    ),
    col_select = gene_id
  )

# List of genes that get activated during transdifferentiation
activated_gene_set <-
  vroom("data/257.notExpressed.0h.txt",
    col_names = "gene_id",
    delim = "\t"
  )

# List of genes that get upregulated during transdifferentiation
# but were already expressed at the beginning
already_expressed_gene_set <-
  vroom("data/629.expressed.0h.txt",
    col_names = "gene_id",
    delim = "\t"
  )

# Gene Classification Table
classes <-
  vroom("data/metadata.class2.tsv",
    skip = 1,
    col_names = c(
      "gene_id",
      "class",
      "time_point",
      "hc",
      "class2",
      "final_class"
    )
  )

# List of genes that get deactivated during transdifferentiation
deactivated_gene_set <-
  classes %>%
  filter(final_class == "downregulation") %>%
  select(gene_id) %>%
  left_join(gene_exp) %>%
  filter(`168` < 1) %>%
  select(gene_id)

# Final Table that includes all classifications
classes <-
  gene_exp %>%
  select(1) %>%
  mutate(
    activated = gene_id %in% activated_gene_set$gene_id,
    already_expressed = gene_id %in% already_expressed_gene_set$gene_id,
    deactivated = gene_id %in% deactivated_gene_set$gene_id,
    stably_expressed =
    !(gene_id %in% de_genes$gene_id) & !(gene_id %in% silenced_genes),
    silent = gene_id %in% silenced_genes
  ) %>%
  left_join(
    classes %>%
      select(gene_id, final_class)
  ) %>%
  pivot_longer(
    activated:silent, names_to = "class", values_to = "class_bool"
  ) %>%
  mutate(
    class = ifelse(!class_bool, final_class, class)
  ) %>%
  drop_na(class) %>%
  select(gene_id, class) %>%
  distinct() %>%
  mutate(class = factor(
    class,
    levels = c(
      "silent", "stably_expressed", "activated", "already_expressed",
      "deactivated", "upregulation", "downregulation", "peaking", "bending"
    ),
    labels = c(
      "Silent", "Stable", "Activated (257)", "Already ON (629)",
      "Deactivated (251)", "Upregulation", "Downregulation", "Peaking", "Bending"
    )
  ))

# Overlap between promoter regions and cCREs
promoter_ccres_overlap <-
  vroom(
    "data/cCREs/promoter_cCREs_intersection.tsv",
    col_names = FALSE,
    col_select = c(6, 7, 11, 14)
  )

# Chromatin state assignment for each cCRE in all timepoints
thmm_profiles <-
  vroom(
    "data/cCREs/HMM.6.gene.matrix.1.tsv",
    col_names = FALSE,
    skip = 1
  )

names(thmm_profiles) <-
  c(
    "cCRE",
    "000", "003", "006", "009", "012", "018",
    "024", "036", "048", "072", "120", "168"
  )


#-------------------------------------------
## Get a unique state profile for each gene
#-------------------------------------------

# Join the `promoter_ccres_overlap` table with the `thmm_profiles` table
# to get a gene assignment for every cCRE with a state profile
thmm_profiles_per_gene <-
  promoter_ccres_overlap %>%
  separate_wider_delim(X14, ".", names = c("gene_id", "version")) %>%
  select(-version) %>%
  inner_join(
    thmm_profiles,
    by = c("X7" = "cCRE")
  )

# Select cCREs with a PLS or pELS classification
# Keep a single state profile per gene - cCRE pair
# independently of the transcript id in `X11` column
thmm_profiles_per_gene_clean <-
  thmm_profiles_per_gene %>%
  filter(X6 %in% c("PLS", "pELS")) %>%
  select(-X11) %>%
  distinct()

# Check the distribution of the number of cCREs per gene
thmm_profiles_per_gene_clean %>%
  group_by(gene_id) %>%
  summarise(n = n()) %>%
  pull(n) %>%
  table()

# The same but checking the number of distinct state profiles instead of the
# number of cCREs
thmm_profiles_per_gene_clean %>%
  select(-c(1:2)) %>%
  distinct() %>%
  group_by(gene_id) %>%
  summarise(n = n()) %>%
  pull(n) %>%
  table()

# Per gene and time-point keep the state assigned to the highest number of
# cCREs associated with a gene (voting scheme)
# In the case of a tie show preference to the bivalent state and then to the
# most "active state"
gene_profiles_longer <-
  thmm_profiles_per_gene_clean %>%
  rename(
    class = X6,
    cCRE_id = X7
  ) %>%
  pivot_longer(
    -c(1:3),
    names_to = "time",
    values_to = "state"
  ) %>%
  mutate(time = as.integer(time)) %>%
  group_by(gene_id, time, state) %>%
  summarise(N = n()) %>%
  filter(N == max(N)) %>%
  summarise(
    state = case_when(
      3 %in% state ~ 3,
      6 %in% state ~ 6,
      5 %in% state ~ 5,
      4 %in% state ~ 4,
      2 %in% state ~ 2,
      1 %in% state ~ 1
    )
  ) %>%
  ungroup()


#------------------------------------
## Check the frequency of transitions
#------------------------------------

# Choose genes whose promoter transition
# to "more active" state without going back
# to a "less active" state at any following
# timepoint (monotonically increasing)
#
# Change "inc" to "dec" in the mon function
# to find transitions to a less active state
# (monotonically decreasing)
n_inc_gene_profiles <-
  gene_profiles_longer %>%
  group_by(gene_id) %>%
  summarise(keep = mon(state, "inc")) %>%
  filter(keep) %>%
  select(-keep) %>%
  left_join(gene_profiles_longer) %>%
  group_by(gene_id) %>%
  summarise(
    n = state %>%
      unique() %>%
      length()
  )

# Check the distribution of the number of states per gene
# across all timepoints during transdifferentiation
n_inc_gene_profiles %>%
  rename(n_states = n) %>%
  group_by(n_states) %>%
  summarise(n_genes = n())

# Control genes are the ones whose promoter remains
# in the state in all timepoints
control_genes <-
  n_inc_gene_profiles %>%
  filter(n == 1) %>%
  select(-n)

# How many controls do we have per state?
control_genes %>%
  left_join(gene_profiles_longer) %>%
  select(gene_id, state) %>%
  distinct() %>%
  group_by(state) %>%
  summarise(n_genes = n())

# How many controls do we have per state and per gene category/class?
control_genes %>%
  left_join(gene_profiles_longer) %>%
  select(gene_id, state) %>%
  distinct() %>%
  left_join(classes) %>%
  group_by(state, class) %>%
  summarise(n_genes = n()) %>%
  print(n = 100)

# Test genes are the ones whose promoter undergoes
# a single transition during transdifferentiation
test_genes <-
  n_inc_gene_profiles %>%
  filter(n == 2) %>%
  select(-n)

# A table that contains a string declaring
# the transition per gene, e.g. 3 -> 6
transitions <-
  test_genes %>%
  left_join(gene_profiles_longer) %>%
  group_by(gene_id) %>%
  summarise(
    transition = state %>%
      unique() %>%
      paste(collapse = " -> ")
  )

# Check the distribution of the number of genes per transition
transitions %>%
  group_by(transition) %>%
  summarise(n_genes = n()) %>%
  arrange(desc(n_genes))

# How many test genes do we have per transition and per gene category/class?
transitions %>%
  left_join(classes) %>%
  group_by(transition, class) %>%
  summarise(n_genes = n()) %>%
  print(n = 100)


#------------------------
## Test 3 -> 6 transition
## (the most frequent)
#------------------------

# A function that fetches data for a specific transition
# That is RNA expression values
# and a binary vector associated to the moment of transition
# It does that for all test genes that undergo a specific transition
# and their respective controls
# e.g. test genes: x -> y, controls remain at state x
get_test_data <- function(from, to){
  control_genes %>%
    left_join(gene_profiles_longer) %>%
    select(gene_id, state) %>%
    distinct() %>%
    filter(state == from) %>%
    select(gene_id) %>%
    bind_rows(
      transitions %>%
        filter(transition == paste(from, to, sep = " -> ")) %>%
        select(gene_id)
    ) %>%
    left_join(gene_profiles_longer) %>%
    left_join(
      gene_exp_requantile %>%
        rename(
          gene_exp = expr_requant,
          time = timepoint
        )
    ) %>%
    mutate(
      state = ifelse(state == from, 0, ifelse(state == to, 1, NA))
    )
}

# Fetch data for 3 -> 6 transition
data_3to6 <-
  get_test_data(3, 6)

# How many controls?
data_3to6 %>%
  select(gene_id) %>%
  distinct() %>%
  inner_join(control_genes)

# How many test genes?
data_3to6 %>%
  select(gene_id) %>%
  distinct() %>%
  inner_join(test_genes)

out_causality_appearance <- list()
trim_npre <- 2

# Run causality using the tjbal function
out_causality_appearance[["3 -> 6"]] <-
  tjbal(
    gene_exp ~ state, data = data_3to6,
    index = c("gene_id", "time"), estimator = "kernel",
    demean = FALSE, vce = "jackknife", nsims = 1000,
    cores = 1, parallel = FALSE, seed = 1234, trim.npre = trim_npre
  )


#-----------------------
## Test more transitions
#-----------------------

# The rest of transitions to test
# The ones for which there is a
# minimum number of test and control genes
to_test <-
  transitions %>%
  group_by(transition) %>%
  summarise(n_genes = n()) %>%
  arrange(desc(n_genes)) %>%
  slice(2:6) %>%
  pull(transition)

# Run causality analysis
for(t in to_test){

  print(t)
  t_split <- strsplit(t, " -> ")
  from <- t_split[[1]][1]
  to <- t_split[[1]][2]

  data_to_test <-
    get_test_data(from, to)

  # control genes
  data_to_test %>%
    select(gene_id) %>%
    distinct() %>%
    inner_join(control_genes) %>%
    print()

  # test genes
  data_to_test %>%
    select(gene_id) %>%
    distinct() %>%
    inner_join(test_genes) %>%
    print()

  out_causality_appearance[[t]] <-
    tjbal(
      gene_exp ~ state, data = data_to_test,
      index = c("gene_id", "time"), estimator = "kernel",
      demean = FALSE, vce = "jackknife", nsims = 1000,
      cores = 1, parallel = FALSE, seed = 1234, trim.npre = trim_npre
    )
}
