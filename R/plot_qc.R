#' Plot QC Summary Dashboard
#'
#' @description Multi-panel QC dashboard combining read counts per bin,
#' classification breakdown, and read statistics using patchwork.
#'
#' @param cys_data Data frame with position and bin columns.
#' @param n_bins Integer. Number of bins.
#' @param bin_labels Character vector. Custom bin labels.
#' @param classification_summary Data frame with mutation class summary
#'   (optional, from \code{\link{classify_mutations}}).
#' @param title Character. Overall dashboard title.
#'
#' @return A patchwork composite ggplot.
#' @export
plot_qc_summary <- function(cys_data,
                              n_bins = 4L,
                              bin_labels = NULL,
                              classification_summary = NULL,
                              title = "Deep Sequencing QC Dashboard") {

  if (is.null(bin_labels)) {
    bin_labels <- paste("Bin", seq_len(n_bins))
  }
  bin_cols <- paste0("bin_", seq_len(n_bins))

  # Panel 1: Total reads per bin (bar chart)
  totals <- data.frame(
    bin = factor(bin_labels, levels = bin_labels),
    total = vapply(bin_cols, function(bc) sum(cys_data[[bc]], na.rm = TRUE),
                   numeric(1)),
    stringsAsFactors = FALSE
  )

  grad_colors <- c("#00d2ff", "#7b2ff7", "#ff6b6b", "#ffd93d",
                   "#6bcb77", "#ee82ee", "#ffa500", "#00ced1")[seq_len(n_bins)]

  p1 <- ggplot2::ggplot(totals, ggplot2::aes(
    x = .data$bin, y = .data$total, fill = .data$bin
  )) +
    ggplot2::geom_col(width = 0.6, alpha = 0.9, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = scales::comma(.data$total)),
                        vjust = -0.5, color = "white", size = 3.5,
                        fontface = "bold") +
    ggplot2::scale_fill_manual(values = grad_colors) +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                 expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = "Total Reads per Bin", x = NULL, y = "Read Count") +
    .dark_theme()

  # Panel 2: Distribution of counts across positions (boxplot per bin)
  long_data <- tidyr::pivot_longer(
    cys_data[, c("position", bin_cols)],
    cols = dplyr::all_of(bin_cols),
    names_to = "bin",
    values_to = "count"
  )
  long_data$bin_label <- factor(
    long_data$bin, levels = bin_cols, labels = bin_labels
  )

  p2 <- ggplot2::ggplot(long_data, ggplot2::aes(
    x = .data$bin_label, y = .data$count + 1, fill = .data$bin_label
  )) +
    ggplot2::geom_boxplot(alpha = 0.7, show.legend = FALSE,
                           outlier.color = "#ffffff60", outlier.size = 0.5) +
    ggplot2::scale_fill_manual(values = grad_colors) +
    ggplot2::scale_y_log10(labels = scales::comma) +
    ggplot2::labs(title = "Read Count Distribution", x = NULL,
                   y = "Count + 1 (log scale)") +
    .dark_theme()

  # Panel 3: Correlation heatmap between bins
  bin_matrix <- as.matrix(cys_data[, bin_cols])
  cor_mat <- cor(bin_matrix, use = "pairwise.complete.obs")
  cor_df <- expand.grid(
    Bin1 = factor(bin_labels, levels = bin_labels),
    Bin2 = factor(bin_labels, levels = rev(bin_labels)),
    stringsAsFactors = FALSE
  )
  cor_df$correlation <- as.vector(cor_mat)

  p3 <- ggplot2::ggplot(cor_df, ggplot2::aes(
    x = .data$Bin1, y = .data$Bin2, fill = .data$correlation
  )) +
    ggplot2::geom_tile(color = "#1a1a2e") +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$correlation)),
      color = "white", size = 3.5, fontface = "bold"
    ) +
    ggplot2::scale_fill_viridis_c(option = "plasma", limits = c(-1, 1),
                                    name = "r") +
    ggplot2::labs(title = "Bin Correlation", x = NULL, y = NULL) +
    .dark_theme() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )

  # Panel 4: Top variable positions
  cys_data$cv <- apply(bin_matrix, 1, function(x) {
    m <- mean(x)
    if (m == 0) return(0)
    sd(x) / m
  })
  top_var <- cys_data[order(-cys_data$cv), ][seq_len(min(15, nrow(cys_data))), ]
  top_var$label <- paste0(top_var$wt_residue, top_var$position)
  top_var$label <- factor(top_var$label, levels = rev(top_var$label))

  p4 <- ggplot2::ggplot(top_var, ggplot2::aes(
    x = .data$label, y = .data$cv
  )) +
    ggplot2::geom_segment(ggplot2::aes(xend = .data$label, y = 0,
                                         yend = .data$cv),
                           color = "#4ecdc4", linewidth = 1) +
    ggplot2::geom_point(color = "#ff6b6b", size = 3) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Most Variable Positions",
                   x = NULL, y = "Coefficient of Variation") +
    .dark_theme()

  # Compose dashboard
  dashboard <- patchwork::wrap_plots(p1, p2, p3, p4, ncol = 2) +
    patchwork::plot_annotation(
      title = title,
      subtitle = "CcdB Cysteine Scanning Mutagenesis Library",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 18,
                                            hjust = 0.5, color = "white"),
        plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey70",
                                               size = 12),
        plot.background = ggplot2::element_rect(fill = "#0f0f23", color = NA)
      )
    )

  dashboard
}


# Internal dark theme helper
.dark_theme <- function() {
  ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12,
                                          hjust = 0.5, color = "white"),
      plot.background = ggplot2::element_rect(fill = "#0f0f23", color = NA),
      panel.background = ggplot2::element_rect(fill = "#0f0f23", color = NA),
      panel.grid.major = ggplot2::element_line(color = "#ffffff10"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "grey80"),
      axis.title = ggplot2::element_text(color = "grey80"),
      legend.text = ggplot2::element_text(color = "grey80"),
      legend.title = ggplot2::element_text(color = "grey80"),
      legend.background = ggplot2::element_rect(fill = "#0f0f23", color = NA)
    )
}
