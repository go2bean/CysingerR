#' Plot Cysteine Enrichment Bar Plot
#'
#' @description Grouped or stacked bar plot of cysteine enrichment per residue
#' position, colored by bin/sort condition.
#'
#' @param cys_data Data frame from \code{\link{aggregate_cys_residue_level}} or
#'   \code{\link{normalize_cys_frequency}}.
#' @param n_bins Integer. Number of bins to display.
#' @param bin_labels Character vector. Custom bin labels.
#' @param normalized Logical. If TRUE, use enrichment columns; if FALSE use
#'   raw bin counts (default FALSE).
#' @param stacked Logical. Stacked bars (TRUE) or grouped (FALSE, default).
#' @param title Character. Plot title.
#' @param highlight_positions Integer vector. Residue positions to highlight
#'   with a different background color.
#'
#' @return A ggplot object.
#' @export
plot_cys_enrichment <- function(cys_data,
                                 n_bins = 4L,
                                 bin_labels = NULL,
                                 normalized = FALSE,
                                 stacked = FALSE,
                                 title = "Cysteine Substitution Read Counts by Position",
                                 highlight_positions = NULL) {

  if (is.null(bin_labels)) {
    bin_labels <- paste("Bin", seq_len(n_bins))
  }

  if (normalized && any(grepl("^enrichment_", names(cys_data)))) {
    value_cols <- grep("^enrichment_", names(cys_data), value = TRUE)
    y_label <- "Enrichment Ratio"
    title <- sub("Read Counts", "Enrichment", title)
  } else {
    value_cols <- paste0("bin_", seq_len(n_bins))
    y_label <- "Read Count"
  }

  plot_data <- cys_data[, c("position", "wt_residue", value_cols)]

  long_data <- tidyr::pivot_longer(
    plot_data,
    cols = dplyr::all_of(value_cols),
    names_to = "bin",
    values_to = "value"
  )

  long_data$bin_label <- factor(
    long_data$bin,
    levels = value_cols,
    labels = if (length(bin_labels) == length(value_cols)) bin_labels
             else paste("Bin", seq_along(value_cols))
  )

  # Premium color palette
  bin_colors <- c("#00d2ff", "#7b2ff7", "#ff6b6b", "#ffd93d",
                  "#6bcb77", "#ee82ee", "#ffa500", "#00ced1")[seq_along(value_cols)]

  position_type <- if (stacked) "stack" else "dodge"

  p <- ggplot2::ggplot(long_data, ggplot2::aes(
    x = .data$position,
    y = .data$value,
    fill = .data$bin_label
  ))

  # Add highlight rectangles if specified
  if (!is.null(highlight_positions)) {
    for (hp in highlight_positions) {
      p <- p + ggplot2::geom_rect(
        xmin = hp - 0.5, xmax = hp + 0.5,
        ymin = -Inf, ymax = Inf,
        fill = "#ffffff10", inherit.aes = FALSE
      )
    }
  }

  p <- p +
    ggplot2::geom_col(position = position_type, width = 0.7, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = bin_colors, name = "Bin") +
    ggplot2::scale_x_continuous(
      breaks = seq(1, max(long_data$position), by = 10),
      expand = ggplot2::expansion(mult = 0.01)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.05))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "CcdB Cysteine Scanning Library",
      x = "Residue Position",
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16,
                                          hjust = 0.5, color = "white"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey70",
                                             size = 11),
      plot.background = ggplot2::element_rect(fill = "#0f0f23", color = NA),
      panel.background = ggplot2::element_rect(fill = "#0f0f23", color = NA),
      panel.grid.major.y = ggplot2::element_line(color = "#ffffff15"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "grey80"),
      axis.title = ggplot2::element_text(color = "grey80", face = "bold"),
      legend.position = "top",
      legend.text = ggplot2::element_text(color = "grey80"),
      legend.title = ggplot2::element_text(color = "grey80", face = "bold"),
      legend.background = ggplot2::element_rect(fill = "#0f0f23", color = NA),
      legend.key = ggplot2::element_rect(fill = "#0f0f23", color = NA),
      plot.margin = ggplot2::margin(15, 15, 10, 15)
    )

  p
}
