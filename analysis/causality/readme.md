# Analysis Description
The code included in the present scripts was used to investigate the causal
relationship between histone post-translational modifications (hPTMs) and RNA
expression during trans-differentiation of the BLaER1 cell line. To do this, we
used [trajectory balancing](https://yiqingxu.org/packages/tjbal/index.html),
a kernel-based reweighting method for causal inference with panel data.

## Data Description
We tested the effect that a chromatin state transition at promoters (or
enhancers) has on RNA expression of their corresponding genes. We also tested
the effect that gene activation has on the deposition of different hPTMs at
different regions of interest (e.g. promoters or enhancers).

Therefore, we used the following data as input for this analysis:
- Normalized RNA expression values at gene level
- Chromatin state assignment of the studied regulatory elements
- Overlap of the studied regulatory elements with promoter (or enhancer) regions
- Normalized ChIP-Seq signal for hPTMs quantified at promoter (or enhancer) regions

## Environment
The environment including all the libraries required for this analysis has been
packaged in a [Docker](https://www.docker.com) image and is available from
[Docker Hub](https://hub.docker.com):

```
docker pull vntasis/r-env:panel_causality-1.0
```

## Scripts Description
- [`causality_analysis_promoters.R`](scripts/causality_analysis_promoters.R):
  Test the causal effect that a chromatin state change at promoters has on RNA
  expression
- [`causality_analysis_reversed_promoters.R`](scripts/causality_analysis_reversed_promoters.R):
  Test the causal effect that gene activation has on the signal of different
  hPTMs at promoters
- [`causality_analysis_enhancers.R`](scripts/causality_analysis_enhancers.R):
  Test the causal effect that a chromatin state change at enhancers has on RNA
  expression
- [`causality_analysis_reversed_enhancers.R`](scripts/causality_analysis_reversed_enhancers.R):
  Test the causal effect that gene activation has on the signal of different
  hPTMs at enhancers
- [`att_line_plots.R`](scripts/att_line_plots.R): Generate ATT (Average
  Treatment effect) and counterfactual lineplots
- [`summary_plot_att_promoters.R`](scripts/summary_plot_att_promoters.R):
  Calculate relative ATT values and Plot Average Causal Effect for chromatin
  marking at promoters
- [`summary_plot_att_enhancers.R`](scripts/summary_plot_att_enhancers.R):
  Calculate relative ATT values and Plot Average Causal Effect for chromatin
  marking at enhancers
- [`theme_publication.R`](scripts/theme_publication.R): Function for theming
  plots generated with `ggplot2`

## Notes
- These scripts were run in an interactive manner. They are not meant to be
  used as executables.
- The trajectory balancing output is also provided in `.RData` files (which are
  also used as input for the `summary_plot_*` scripts.
- The `causality_analysis_*` scripts do not cover all scenarios. Check the
  comments for cases for which changing a parameter or an input is appropriate.
  For instance, in `causality_analysis_promoters.R`, changing "inc" to "dec" in
  the `mon` function will switch from detecting transitions to a more active
  state, to finding transitions to a less active state.
- Final figures were mounted using [Inkscape](https://inkscape.app).
