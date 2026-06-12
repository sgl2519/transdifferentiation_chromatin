#!/usr/bin/env r

#*******************************************
## Generate ATT and counterfactual lineplots
#*******************************************

#---------------------
## Plotting parameters
#---------------------

font_small <- 0.66
font_med <- 0.69
font_big <- 0.55


#----------------
## Generate plots
#----------------

# Save in a pdf file
# Here is an example from the promoter tests
pdf(
  "plots/causality_transition_to_active_state_promoters.pdf",
  family = "NimbusSan",
  width = 2,
  height = 3
)

# out_causality contains a list whose elements
# are an output of the tjbal function
for (m in names(out_causality)) {
  print(m)

  # Show the number of controls and test genes
  # involved in a test
  plot.new()
  text(
    .5, .5,
    paste0(
      c(
        m,
        paste0("Nco: ", out_causality[[m]]$Nco),
        paste0("Ntr: ", out_causality[[m]]$Ntr)
      ),
      collapse = "\n"
    )
  )

  # Choose the relative timepoints to show in x-axis
  # as the ones that involve at least 10 test genes
  tps <-
    (
     out_causality[[m]]$sub.ntr %>%
       rownames
    )[out_causality[[m]]$sub.ntr %>% rowSums %>% `>=`(10)]

  xlims <-
    c(tps[1], tps[length(tps)]) %>%
    as.integer

  # make plots
  # ATT plot
  plot(
    out_causality[[m]], cex.main = font_big, cex.axis = font_small,
    cex.text = font_small, cex.lab = font_med, cex.legend = font_med, xlim = xlims
  )

  # Counterfactual lineplot
  plot(
    out_causality[[m]], type = "ct", cex.main = font_big,
    cex.axis = font_small, cex.text = font_small, cex.lab = font_med,
    cex.legend = font_med, legend.ncol = 1, xlim = xlims
  )

}

dev.off()
