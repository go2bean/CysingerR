#' Classify Mutations from Substitution Table
#'
#' @description Categorizes reads as WT, single, double, triple, or
#' multi-mutant based on the number of unique residue positions with
#' substitutions. Replaces \code{classify_substitutions.pl}.
#'
#' @param sub_table Data frame from \code{\link{build_substitution_table}}.
#'
#' @return A list with:
#'   \describe{
#'     \item{wt}{Data frame of WT reads.}
#'     \item{single}{Data frame of single-mutant reads.}
#'     \item{multi}{Data frame of multi-mutant reads (double, triple, >3).}
#'     \item{summary}{Data frame with counts for each mutation class.}
#'     \item{wt_per_bin}{Data frame of WT read counts per bin.}
#'   }
#' @export
classify_mutations <- function(sub_table) {
  n <- nrow(sub_table)
  category <- character(n)

  for (i in seq_len(n)) {
    if (sub_table$is_wt[i]) {
      category[i] <- "wt"
    } else {
      # Count unique residue positions
      res_pos <- sub_table$residue_positions[i]
      if (res_pos == "wt") {
        category[i] <- "wt"
      } else {
        unique_pos <- unique(as.integer(strsplit(res_pos, "\\s+")[[1]]))
        num_unique <- length(unique_pos)
        category[i] <- switch(as.character(num_unique),
                              "1" = "single",
                              "2" = "double",
                              "3" = "triple",
                              "multi")
      }
    }
  }

  sub_table$mutation_class <- category

  wt_df <- sub_table[category == "wt", ]
  single_df <- sub_table[category == "single", ]
  multi_df <- sub_table[category %in% c("double", "triple", "multi"), ]

  # Summary counts
  counts <- table(factor(category,
                         levels = c("wt", "single", "double", "triple", "multi")))
  summary_df <- data.frame(
    class = names(counts),
    count = as.integer(counts),
    fraction = as.numeric(counts) / sum(counts),
    stringsAsFactors = FALSE
  )

  # WT counts per bin
  wt_per_bin <- as.data.frame(table(wt_df$bin), stringsAsFactors = FALSE)
  names(wt_per_bin) <- c("bin", "wt_count")

  list(
    wt = wt_df,
    single = single_df,
    multi = multi_df,
    all = sub_table,
    summary = summary_df,
    wt_per_bin = wt_per_bin
  )
}
