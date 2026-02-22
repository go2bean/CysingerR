#' Plot Cysteine Frequency Heatmap
#'
#' @description Publication-quality heatmap of cysteine substitution
#' frequencies across residue positions and bins. Uses ggplot2 with
#' viridis color scale.
#'
#' @param cys_data Data frame from \code{\link{aggregate_cys_residue_level}} or
#'   \code{\link{load_example_data}}.
#' @param n_bins Integer. Number of bins to display.
#' @param bin_labels Character vector. Custom bin labels (e.g., c("Input",
#'   "Unsorted", "Low Sort", "High Sort")).
#' @param log_transform Logical. Apply log10(x+1) transformation (default TRUE).
#' @param title Character. Plot title.
#' @param color_palette Character. Viridis palette name (default "inferno").
#'
#' @return A ggplot object.
#' @export
#' @examples
#' \dontrun{
#' cys_data <- load_example_data("residue_level")
#' plot_cys_heatmap(cys_data)
#' }
plot_cys_heatmap <- function(cys_data,
                              n_bins = 4L,
                              bin_labels = NULL,
                              log_transform = TRUE,
                              title = "Cysteine Substitution Frequency Heatmap",
                              color_palette = "inferno") {
  bin_cols <- paste0("bin_", seq_len(n_bins))

  if (is.null(bin_labels)) {
    bin_labels <- paste("Bin", seq_len(n_bins))
  }

  # Pivot to long format
  plot_data <- cys_data[, c("position", "wt_residue", bin_cols)]
  long_data <- tidyr::pivot_longer(
    plot_data,
    cols = dplyr::all_of(bin_cols),
    names_to = "bin",
    values_to = "count"
  )
  long_data$bin_label <- factor(
    long_data$bin,
    levels = bin_cols,
    labels = bin_labels
  )

  if (log_transform) {
    long_data$value <- log10(long_data$count + 1)
    fill_label <- expression(log[10](count + 1))
  } else {
    long_data$value <- long_data$count
    fill_label <- "Read Count"
  }

  # Create residue labels
  long_data$res_label <- paste0(long_data$wt_residue, long_data$position)

  p <- ggplot2::ggplot(long_data, ggplot2::aes(
    x = .data$position,
    y = .data$bin_label,
    fill = .data$value
  )) +
    ggplot2::geom_tile(color = "grey20", linewidth = 0.1) +
    ggplot2::scale_fill_viridis_c(
      option = color_palette,
      name = fill_label,
      guide = ggplot2::guide_colorbar(
        barwidth = 15, barheight = 0.8,
        title.position = "top", title.hjust = 0.5
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(1, max(long_data$position), by = 10),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "CcdB Cysteine Scanning Mutagenesis Library",
      x = "Residue Position",
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16,
                                          hjust = 0.5, color = "white"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey70",
                                             size = 11),
      plot.background = ggplot2::element_rect(fill = "#1a1a2e", color = NA),
      panel.background = ggplot2::element_rect(fill = "#1a1a2e", color = NA),
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "grey80"),
      axis.title = ggplot2::element_text(color = "grey80", face = "bold"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(color = "grey80"),
      legend.title = ggplot2::element_text(color = "grey80"),
      legend.background = ggplot2::element_rect(fill = "#1a1a2e", color = NA),
      plot.margin = ggplot2::margin(15, 15, 15, 15)
    )

  p
}
