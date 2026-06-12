#!/usr/bin/env r

#******************************************************************
## Calculate relative ATT values and
## Plot Average Causal Effect (both in absolute and relative terms)
## Enahncers version
#******************************************************************

#-------------------------
## Load required libraries
#-------------------------
library(tidyverse)


#-----------------------------------
## Load output of causality analysis
#-----------------------------------
load("RData/final/tjbal_out_appearance_states_enhancers_proximal.RData")
load("RData/final/tjbal_out_disappearance_states_enhancers_proximal.RData")

proximal_appearance_states2exp <- out_causality_appearance
proximal_disappearance_states2exp <- out_causality_disappearance

load("RData/final/tjbal_out_appearance_states_enhancers_rE2G.RData")
load("RData/final/tjbal_out_disappearance_states_enhancers_rE2G.RData")

rE2G_appearance_states2exp <- out_causality_appearance
rE2G_disappearance_states2exp <- out_causality_disappearance

load("RData/final/tjbal_out_appearance_gene_exp_as_treatment_enhancers_proximal.RData")

proximal_appearance_exp2marks <- out_causality_appearance

load("RData/final/tjbal_out_appearance_gene_exp_as_treatment_enhancers_rE2G.RData")

rE2G_appearance_exp2marks <- out_causality_appearance

rm(list = c("out_causality_appearance", "out_causality_disappearance"))

#----------------
## Get ATT values
#----------------
ATT <- tibble()

for(data in c(
  "proximal_appearance_states2exp", "proximal_disappearance_states2exp",
  "rE2G_appearance_states2exp", "rE2G_disappearance_states2exp",
  "proximal_appearance_exp2marks", "rE2G_appearance_exp2marks"
)) {

  assign("causality_data", get(data))
  label <- data

  # ATT and P-values for all tests
  att <-
    causality_data %>%
    map(
      \(x) as.data.frame(x$est.att.avg) %>% as_tibble()
    ) %>%
    list_rbind(names_to = "treatment") %>%
    select(treatment, ATT, p.value) %>%
    rename(att = ATT)

  # Number of treated and control genes
  Ntr_Nco <-
    causality_data %>%
    map(
      \(x) tibble(Ntr = x$Ntr, Nco = x$Nco)
    ) %>%
    list_rbind(names_to = "treatment")

  # Merge all data
    ATT <-
      att %>%
      full_join(Ntr_Nco) %>%
      mutate(test = label) %>%
      bind_rows(ATT)

}


# Adjust p-values for multiple testing
ATT <-
  ATT %>%
  mutate(sign = ifelse(att >= 0, "Positive Effect", "Negative Effect")) %>%
  group_by(test) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH")) %>%
  ungroup()


#------------------------
## Calculate relative ATT
#------------------------
rel_ATT <- tibble()

for(data in c(
  "proximal_appearance_states2exp", "proximal_disappearance_states2exp",
  "rE2G_appearance_states2exp", "rE2G_disappearance_states2exp",
  "proximal_appearance_exp2marks", "rE2G_appearance_exp2marks"
)) {

  assign("causality_data", get(data))
  label <- data

  # Average response values for all relative timepoints and tests
  Y_bar <-
    causality_data %>%
    map(
      \(x) as.data.frame(x$Y.bar) %>% as_tibble(rownames = "timepoint")
    ) %>%
    list_rbind(names_to = "treatment")

  # Number of treated genes per relative timepoint
  N_treated <-
    causality_data %>%
    map(
      \(x) as.data.frame(x$ntreated) %>% as_tibble(rownames = "timepoint")
    ) %>%
    list_rbind(names_to = "treatment") %>%
    rename(N = 3)


  # Merge response variable data with number of treated cases
  # Categorise data-points to pre- and post-treatment
  Y_N <-
    Y_bar %>%
    full_join(N_treated) %>%
    mutate(
      timepoint = as.integer(timepoint),
      period = ifelse(timepoint <= 0, "pre", "post")
    )


  # Calculate relative ATT
  Y_summary <-
    Y_N %>%
    group_by(treatment, period) %>%
    summarise(
      Y_tr = mean(Y.bar.tr),
      Y_ct = mean(Y.bar.ct),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = period, values_from = c(Y_tr, Y_ct)) %>%
    mutate(att_rel = log2((Y_tr_post / Y_tr_pre) / (Y_ct_post / Y_ct_pre)))

  # Calculate relative ATT (weighted version)
  Y_weighted_summary <-
    Y_N %>%
    group_by(treatment, period) %>%
    summarise(
      Y_tr = weighted.mean(x = Y.bar.tr, w = N),
      Y_ct = weighted.mean(x = Y.bar.ct, w = N),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = period, values_from = c(Y_tr, Y_ct)) %>%
    mutate(weighted_att_rel = log2((Y_tr_post / Y_tr_pre) / (Y_ct_post / Y_ct_pre)))

  # Final output matrix
  final <-
    Y_summary %>%
    select(treatment, att_rel) %>%
    full_join(
      Y_weighted_summary %>%
        select(treatment, weighted_att_rel)
    ) %>%
    mutate(test = label)

    rel_ATT <- bind_rows(rel_ATT, final)

}



#--------------------------
## Merge the two ATT tables
#--------------------------
ATT_final <-
  rel_ATT %>%
  full_join(ATT) %>%
  select(test, treatment, Ntr, Nco, everything())



#--------------------
## Plotting Variables
#--------------------
source("scripts/theme_publication.R")

treatments <-
  c(
    "1 -> 2", "2 -> 4",
    "4 -> 5", "5 -> 6",
    "3 -> 5", "3 -> 6",

    "6 -> 5", "5 -> 4",
    "4 -> 2", "2 -> 1",
    "6 -> 3", "5 -> 3",

    "H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3",
    "H3K36me3", "H4K20me1", "H3K27me3", "H3K9me3"
  )

treatment_labels <-
  c(
    "Absent to Basal", "Basal to Intermediate",
    "Intermediate to Prominent", "Prominent to Strong",
    "Bivalent to Prominent", "Bivalent to Strong",

    "Strong to Prominent", "Prominent to Intermediate",
    "Intermediate to Basal", "Basal to Absent",
    "Strong to Bivalent", "Prominent to Bivalent",

    "H3K4me1", "H3K4me2", "H3K27ac", "H3K9ac", "H3K4me3",
    "H3K36me3", "H4K20me1", "H3K27me3", "H3K9me3"
  )

treatment_colors <-
  c(
    "#a6cee3", "#1f78b4",
    "#b2df8a", "#33a02c",
    "#fb9a99", "#e31a1c",

    "#fdbf6f", "#ff7f00",
    "#cab2d6", "#6a3d9a",
    "#ffff99", "#b15928",

    "#E5AB00", "#A67C00", "#630039", "#AF4D85", "#D199B9",
    "#7FBC41", "#4C7027", "#1D2976", "#A7ADD4"
  )

treatment_colors_named <- treatment_colors
names(treatment_colors_named) <- treatment_labels

tests <-
  c(
  "proximal_appearance_states2exp", "proximal_disappearance_states2exp",
  "rE2G_appearance_states2exp", "rE2G_disappearance_states2exp",
  "proximal_appearance_exp2marks", "rE2G_appearance_exp2marks"
  )

test_labels <-
  c(
    "Treatment: State Transition in proximal enhancers\nOutcome: Gene Expression\nTransition to more active state\n",
    "Treatment: State Transition in proximal enhancers\nOutcome: Gene Expression\nTransition to less active state\n",
    "Treatment: State Transition in distal enhancers\nOutcome: Gene Expression\nTransition to more active state\n",
    "Treatment: State Transition in distal enhancers\nOutcome: Gene Expression\nTransition to less active state\n",
    "Treatment: Gene Activation\nOutcome: Histone PTM signal in proximal enhancers\n",
    "Treatment: Gene Activation\nOutcome: Histone PTM signal in distal enhancers\n"
  )

measures <- c("att", "weighted_att_rel")
measure_labels <- c("Measure: ATT", "Measure: Relative ATT")


#-----------
## Make Plot
#-----------

# Reshape ATT table
# For state transition consider all cases present also at the promoter tests
# Use factors to fix the order and adjust labels
# Create significance level variable
# Create group variable for splitting the data
to_plot <-
  tibble(treatment = treatments, include = 1) %>%
  head(6) %>%
  left_join(
    tibble(test = tests, include = 1) %>%
      slice(1, 3),
    relationship = "many-to-many"
  ) %>%
  select(-include) %>%
  add_row(
    tibble(treatment = treatments, include = 1) %>%
    slice(7:12) %>%
    left_join(
      tibble(test = tests, include = 1) %>%
        slice(2, 4),
      relationship = "many-to-many"
    ) %>%
    select(-include)
  ) %>%
  full_join(ATT_final) %>%
  select(-c("att_rel", "p.value", "sign")) %>%
  replace_na(list(att = 0, weighted_att_rel = 0)) %>%
  pivot_longer(
    c(att, weighted_att_rel),
    names_to = "measure",
    values_to = "effect"
  ) %>%
  mutate(
    group = str_extract(test, "[^_]+_[^_]+$"),
    test = factor(test, levels = tests, labels = test_labels),
    treatment = factor(treatment, levels = treatments, labels = treatment_labels),
    measure = factor(measure, levels = measures, labels = measure_labels),
    sig.level = case_when(
      is.na(p.adj) ~ "N/A",
      p.adj <= 0.001 ~ "***",
      p.adj <= 0.01 ~ "**",
      p.adj <= 0.05 ~ "*",
      .default = "ns"
    )
  ) %>%
  rename(
    Ntreated = Ntr,
    Ncontrol = Nco
  )

# Split data to make pairs of plot
# with proximal and distal enhancer sharing y-axis
to_plot_split <-
  to_plot %>%
  group_by(group, measure) %>%
  group_split() %>%
  map(
    \(x) x %>%
      group_by(treatment, test, measure) %>%
      group_trim() %>%
      ungroup %>%
      select(-group)
  )

# ATT barplot - split version
pdf(
  "plots/ATT_summary_enhancers_split.pdf",
  width = 2.3, height = 4
)

to_plot_split[-c(5:6)] %>%
  walk(
    function(x) {
      print(
        x %>%
          ggplot(aes(x = treatment, y = effect, fill = treatment)) +
          geom_bar(stat = "identity", color = "black") +
          geom_text(aes(label = sig.level), vjust = "outward", size = 3) +
          facet_wrap(
            ~ test + measure,
            ncol = 1,
            labeller = labeller(
              test = "label_value",
              measure = "label_value",
              .multi_line = FALSE
            )
          ) +
          theme_publication(
            font = "Helvetica",
            size_medium = 10,
            size_small = 9
          ) +
          labs(
            y = "Estimated Causal Effect",
            x = "",
            fill = "",
          ) +
          theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = "right"
          ) +
          guides(fill = "none") +
          scale_fill_manual(
            values = x %>%
              select(treatment) %>%
              arrange(treatment) %>%
              distinct() %>%
              pull(treatment) %>%
              as.character() %>%
              `[`(treatment_colors_named, .) %>%
              unname()
          ) +
          coord_cartesian(clip = "off") +
          scale_y_continuous(
            labels = scales::number_format(accuracy = 0.1),
            breaks = scales::breaks_width(
              (function(df) {
                m <- df %>%
                  pull(effect) %>%
                  abs() %>%
                  max()

                ifelse(m > 2, 1,
                  ifelse(m > 1, 0.5,
                    ifelse(m > 0.5, 0.2, 0.1)
                  )
                )
              }) (x)
            )
          )
        )
    }
  )

dev.off()

# ATT barplot
pdf(
  "plots/ATT_summary_enhancers.pdf",
  width = 12, height = 10
)

to_plot %>%
  ggplot(aes(x = treatment, y = effect, fill = treatment)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = sig.level, vjust = effect > 0), size = 3) +
  facet_wrap(
    ~ test + measure,
    scales = "free", ncol = 3,
    labeller = labeller(
      test = "label_value",
      measure = "label_value",
      .multi_line = FALSE
    )
  ) +
  theme_publication(
    font = "Helvetica",
    size_medium = 10,
    size_small = 9
  ) +
  labs(
    y = "Estimated Causal Effect",
    x = "",
    fill = "",
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right"
  ) +
  guides(fill = guide_legend(nrow = 12)) +
  scale_fill_manual(values = treatment_colors) +
  scale_y_continuous(breaks = scales::breaks_extended(4)) +
  coord_cartesian(clip = "off")

dev.off()
