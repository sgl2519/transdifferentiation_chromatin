theme_publication <- function(
  size_small = 8, size_medium = 9, size_big = 10,
  frame_size = 1, font = "NimbusSan"
) {
  (theme_bw() +
    theme(
      panel.border = element_rect(linewidth = frame_size),
      text = element_text(size = size_big, family = font),
      axis.text = element_text(size = size_small),
      legend.text = element_text(size = size_medium),
      strip.text = element_text(size = size_big, face = "bold"),
      strip.background = element_blank(),
      plot.tag = element_text(face = "bold")
    )
  )
}
