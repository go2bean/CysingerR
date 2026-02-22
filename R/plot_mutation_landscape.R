#' Plot Mutation Landscape
#'
#' @description Lollipop/Manhattan-style plot showing mutation frequency along
#' the protein sequence.
#'
#' @param cys_data Data frame with at least \code{position} and bin count
#'   columns.
#' @param n_bins Integer. Number of bins.
#' @param bin_to_plot Integer. Which bin column to plot (default: highest bin).
#' @param bin_label Character. Label for the plotted bin.
#' @param threshold Numeric. Highlight positions above this count threshold.
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#' @export
plot_mutation_landscape <- function(cys_data,
                                     n_bins = 4L,
                                     bin_to_plot = NULL,
                                     bin_label = "Highest Sort",
                                     threshold = NULL,
                                     title = "Cysteine Mutation Landscape") {
  if (is.null(bin_to_plot)) bin_to_plot <- n_bins
  value_col <- paste0("bin_", bin_to_plot)

  plot_data <- data.frame(
    position = cys_data$position,
    wt_residue = cys_data$wt_residue,
    count = cys_data[[value_col]],
    stringsAsFactors = FALSE
  )

  # Assign highlight category
  if (is.null(threshold)) {
    threshold <- quantile(plot_data$count, 0.75)
  }
  plot_data$highlight <- ifelse(plot_data$count >= threshold, "High", "Low")

  colors <- c("High" = "#ff6b6b", "Low" = "#4ecdc4")

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$position, y = .data$count
  )) +
    ggplot2::geom_segment(
      ggplot2::aes(xend = .data$position, y = 0, yend = .data$count,
                    color = .data$highlight),
      linewidth = 0.6, alpha = 0.7
    ) +
    ggplot2::geom_point(
      ggplot2::aes(color = .data$highlight, size = .data$count),
      alpha = 0.85
    ) +
    ggplot2::scale_color_manual(values = colors, name = "Category",
                                 guide = "none") +
    ggplot2::scale_size_continuous(range = c(0.8, 4), guide = "none") +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed",
                         color = "#ffd93d", alpha = 0.6) +
    ggplot2::annotate("text", x = max(plot_data$position) * 0.95,
                       y = threshold * 1.15,
                       label = paste0("Threshold: ", round(threshold)),
                       color = "#ffd93d", size = 3.5, hjust = 1) +
    # Label top hits
    ggplot2::geom_text(
      data = plot_data[plot_data$count >= quantile(plot_data$count, 0.9), ],
      ggplot2::aes(label = paste0(.data$wt_residue, .data$position)),
      vjust = -1.2, size = 2.8, color = "#ffffff90", fontface = "bold"
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(1, max(plot_data$position), by = 10),
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::scale_y_continuous(
      labels = scales::comma,
      expand = ggplot2::expansion(mult = c(0, 0.12))
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste("Bin:", bin_label),
      x = "Residue Position",
      y = "Read Count"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16,
                                          hjust = 0.5, color = "white"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#4ecdc4",
                                             size = 11, face = "italic"),
      plot.background = ggplot2::element_rect(fill = "#16213e", color = NA),
      panel.background = ggplot2::element_rect(fill = "#16213e", color = NA),
      panel.grid.major.y = ggplot2::element_line(color = "#ffffff10"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "grey80"),
      axis.title = ggplot2::element_text(color = "grey80", face = "bold"),
      plot.margin = ggplot2::margin(15, 20, 10, 15)
    )

  p
}


#' Plot Mutation Class Donut Chart
#'
#' @description Elegant donut chart showing the proportions of WT, single,
#' double, triple, and multi-mutant reads.
#'
#' @param classification List from \code{\link{classify_mutations}} (uses
#'   the \code{summary} element), OR a data frame with columns \code{class}
#'   and \code{count}.
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#' @export
plot_mutation_class_donut <- function(classification,
                                       title = "Mutation Classification") {
  if (is.list(classification) && !is.data.frame(classification)) {
    summary_df <- classification$summary
  } else {
    summary_df <- classification
  }

  summary_df$fraction <- summary_df$count / sum(summary_df$count)
  summary_df$ymax <- cumsum(summary_df$fraction)
  summary_df$ymin <- c(0, summary_df$ymax[-nrow(summary_df)])
  summary_df$ymid <- (summary_df$ymin + summary_df$ymax) / 2
  summary_df$label <- paste0(
    summary_df$class, "\n",
    scales::comma(summary_df$count), "\n",
    sprintf("%.1f%%", summary_df$fraction * 100)
  )

  class_colors <- c(
    "wt" = "#4ecdc4", "single" = "#ff6b6b",
    "double" = "#ffd93d", "triple" = "#7b2ff7",
    "multi" = "#ee82ee"
  )

  p <- ggplot2::ggplot(summary_df, ggplot2::aes(
    ymax = .data$ymax, ymin = .data$ymin,
    xmax = 4, xmin = 2.5,
    fill = .data$class
  )) +
    ggplot2::geom_rect(color = "#0f0f23", linewidth = 0.5) +
    ggplot2::geom_text(
      ggplot2::aes(x = 3.25, y = .data$ymid, label = .data$label),
      color = "white", size = 3, fontface = "bold", lineheight = 0.9
    ) +
    ggplot2::coord_polar(theta = "y") +
    ggplot2::xlim(c(1, 4.5)) +
    ggplot2::scale_fill_manual(values = class_colors, guide = "none") +
    ggplot2::annotate("text", x = 1, y = 0, label = "Mutations",
                       size = 5, color = "white", fontface = "bold") +
    ggplot2::labs(title = title) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16,
                                          hjust = 0.5, color = "white"),
      plot.background = ggplot2::element_rect(fill = "#0f0f23", color = NA),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )

  p
}
