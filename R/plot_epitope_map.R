#' Plot Epitope Profile
#'
#' @description Line plot of normalized enrichment ratio across protein
#' sequence positions — the key figure for mapping binding sites and
#' conformational epitopes, as published in Najar et al. (2017).
#'
#' @param norm_data Data frame from \code{\link{normalize_cys_frequency}}.
#' @param enrichment_bins Integer vector. Which sorted bins to plot (default
#'   all available).
#' @param bin_labels Character vector. Custom labels.
#' @param show_log2 Logical. Plot log2 enrichment instead of raw enrichment
#'   (default TRUE).
#' @param smooth Logical. Add a loess smoother (default TRUE).
#' @param span Numeric. Loess span parameter (default 0.15).
#' @param threshold_line Numeric. Y-value for a significance threshold line.
#' @param highlight_regions List of named vectors with \code{start} and
#'   \code{end} for secondary structure or binding site regions to highlight.
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#' @export
plot_epitope_profile <- function(norm_data,
                                  enrichment_bins = NULL,
                                  bin_labels = NULL,
                                  show_log2 = TRUE,
                                  smooth = TRUE,
                                  span = 0.15,
                                  threshold_line = NULL,
                                  highlight_regions = NULL,
                                  title = "Conformational Epitope Mapping Profile") {

  if (show_log2) {
    enrich_cols <- grep("^log2_enrichment_", names(norm_data), value = TRUE)
    y_label <- expression(log[2]~"Enrichment Ratio")
  } else {
    enrich_cols <- grep("^enrichment_bin_", names(norm_data), value = TRUE)
    y_label <- "Enrichment Ratio"
  }

  if (!is.null(enrichment_bins)) {
    if (show_log2) {
      enrich_cols <- paste0("log2_enrichment_bin_", enrichment_bins)
    } else {
      enrich_cols <- paste0("enrichment_bin_", enrichment_bins)
    }
    enrich_cols <- enrich_cols[enrich_cols %in% names(norm_data)]
  }

  if (is.null(bin_labels)) {
    bin_labels <- paste("Sorted Bin", seq_along(enrich_cols))
  }

  plot_data <- norm_data[, c("position", "wt_residue", enrich_cols)]

  long_data <- tidyr::pivot_longer(
    plot_data,
    cols = dplyr::all_of(enrich_cols),
    names_to = "bin",
    values_to = "enrichment"
  )
  long_data$bin_label <- factor(
    long_data$bin, levels = enrich_cols,
    labels = if (length(bin_labels) == length(enrich_cols)) bin_labels
             else paste("Bin", seq_along(enrich_cols))
  )

  line_colors <- c("#00d2ff", "#ff6b6b", "#ffd93d", "#7b2ff7",
                   "#6bcb77", "#ee82ee")[seq_along(enrich_cols)]

  p <- ggplot2::ggplot(long_data, ggplot2::aes(
    x = .data$position, y = .data$enrichment,
    color = .data$bin_label, group = .data$bin_label
  ))

  # Add highlight regions
  if (!is.null(highlight_regions)) {
    for (region in highlight_regions) {
      label <- if (!is.null(region$label)) region$label else ""
      p <- p + ggplot2::annotate(
        "rect",
        xmin = region$start, xmax = region$end,
        ymin = -Inf, ymax = Inf,
        fill = "#ffffff08", alpha = 0.3
      )
      if (nchar(label) > 0) {
        p <- p + ggplot2::annotate(
          "text",
          x = (region$start + region$end) / 2,
          y = max(long_data$enrichment, na.rm = TRUE) * 0.95,
          label = label, color = "#ffffff60", size = 2.8, angle = 90
        )
      }
    }
  }

  # Add threshold line
  if (!is.null(threshold_line)) {
    p <- p + ggplot2::geom_hline(
      yintercept = threshold_line,
      linetype = "dashed", color = "#ffd93d", alpha = 0.5
    )
  }

  p <- p +
    ggplot2::geom_line(alpha = 0.3, linewidth = 0.4) +
    ggplot2::geom_point(size = 1.2, alpha = 0.5)

  if (smooth) {
    p <- p + ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = TRUE,
      span = span, alpha = 0.15, linewidth = 1.2
    )
  }

  p <- p +
    ggplot2::scale_color_manual(values = line_colors, name = "Sort Condition") +
    ggplot2::scale_x_continuous(
      breaks = seq(1, max(long_data$position), by = 10),
      expand = ggplot2::expansion(mult = 0.02)
    ) +
    ggplot2::labs(
      title = title,
      subtitle = "Cysteine Labeling + Deep Sequencing",
      x = "Residue Position",
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16,
                                          hjust = 0.5, color = "white"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#4ecdc4",
                                             face = "italic", size = 11),
      plot.background = ggplot2::element_rect(fill = "#0a0a1a", color = NA),
      panel.background = ggplot2::element_rect(fill = "#0a0a1a", color = NA),
      panel.grid.major = ggplot2::element_line(color = "#ffffff08"),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(color = "grey80"),
      axis.title = ggplot2::element_text(color = "grey80", face = "bold"),
      legend.position = "top",
      legend.text = ggplot2::element_text(color = "grey80"),
      legend.title = ggplot2::element_text(color = "grey80", face = "bold"),
      legend.background = ggplot2::element_rect(fill = "#0a0a1a", color = NA),
      legend.key = ggplot2::element_rect(fill = "#0a0a1a", color = NA),
      plot.margin = ggplot2::margin(15, 20, 10, 15)
    )

  p
}


#' Plot Surface Accessibility Radar
#'
#' @description Polar/radar-style plot showing cysteine labeling accessibility
#' around the protein sequence. Positions with high enrichment indicate
#' surface-accessible sites.
#'
#' @param norm_data Data frame from \code{\link{normalize_cys_frequency}}.
#' @param enrichment_bin Integer. Which sorted bin enrichment to use.
#' @param show_log2 Logical. Use log2 enrichment (default TRUE).
#' @param title Character. Plot title.
#'
#' @return A ggplot object.
#' @export
plot_surface_accessibility <- function(norm_data,
                                        enrichment_bin = 4L,
                                        show_log2 = TRUE,
                                        title = "Surface Accessibility Map") {
  if (show_log2) {
    col_name <- paste0("log2_enrichment_bin_", enrichment_bin)
  } else {
    col_name <- paste0("enrichment_bin_", enrichment_bin)
  }

  plot_data <- data.frame(
    position = norm_data$position,
    wt_residue = norm_data$wt_residue,
    value = norm_data[[col_name]],
    stringsAsFactors = FALSE
  )

  # Classify by enrichment
  med_val <- median(plot_data$value, na.rm = TRUE)
  plot_data$category <- ifelse(plot_data$value >= med_val, "Accessible", "Buried")

  colors <- c("Accessible" = "#ff6b6b", "Buried" = "#4ecdc4")

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$position, y = .data$value, fill = .data$category
  )) +
    ggplot2::geom_col(width = 1, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = colors, name = "Classification") +
    ggplot2::coord_polar(start = 0) +
    ggplot2::labs(
      title = title,
      subtitle = paste("Bin", enrichment_bin),
      x = "Residue Position",
      y = if (show_log2) expression(log[2]~"Enrichment") else "Enrichment"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 16,
                                          hjust = 0.5, color = "white"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "#4ecdc4",
                                             size = 10),
      plot.background = ggplot2::element_rect(fill = "#0a0a1a", color = NA),
      panel.background = ggplot2::element_rect(fill = "#0a0a1a", color = NA),
      panel.grid = ggplot2::element_line(color = "#ffffff10"),
      axis.text = ggplot2::element_text(color = "grey60", size = 7),
      axis.title = ggplot2::element_text(color = "grey70"),
      legend.position = "bottom",
      legend.text = ggplot2::element_text(color = "grey80"),
      legend.title = ggplot2::element_text(color = "grey80", face = "bold"),
      legend.background = ggplot2::element_rect(fill = "#0a0a1a", color = NA),
      legend.key = ggplot2::element_rect(fill = "#0a0a1a", color = NA),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )

  p
}
