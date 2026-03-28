#' Default lineagefreq color palette (colorblind-safe, Okabe-Ito)
#' @noRd
.lfq_palette <- function(n) {
  base <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999",
    "#000000"
  )
  if (n <= length(base)) base[seq_len(n)]
  else grDevices::colorRampPalette(base)(n)
}

#' Fill scale ?"Other" always gray
#' @noRd
.lfq_scale_fill <- function(variants) {
  cols <- .lfq_palette(length(variants))
  names(cols) <- variants
  if ("Other" %in% variants) cols["Other"] <- "#BBBBBB"
  ggplot2::scale_fill_manual(values = cols)
}

#' Colour scale ?"Other" always gray
#' @noRd
.lfq_scale_colour <- function(variants) {
  cols <- .lfq_palette(length(variants))
  names(cols) <- variants
  if ("Other" %in% variants) cols["Other"] <- "#BBBBBB"
  ggplot2::scale_colour_manual(values = cols)
}

#' Publication-ready theme
#' @noRd
.lfq_theme <- function(base_size = 12) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.title     = ggplot2::element_text(face = "bold",
                                               size = base_size - 1),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title       = ggplot2::element_text(size = base_size),
      strip.text       = ggplot2::element_text(face = "bold")
    )
}
