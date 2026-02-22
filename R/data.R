#' Load Example Cysteine Frequency Data
#'
#' @description Loads the built-in CcdB cysteine frequency data from the
#' sample output of the original pipeline.
#'
#' @param type Character. One of \code{"residue_level"} (default) to load
#'   residue-level cysteine frequencies, or \code{"wt_sequence"} to load
#'   the wild-type FASTA.
#'
#' @return Data frame (for frequency data) or character string (for WT
#'   sequence).
#' @export
#' @examples
#' # Load the sample residue-level data
#' cys_data <- load_example_data("residue_level")
#' head(cys_data)
load_example_data <- function(type = c("residue_level", "wt_sequence")) {
  type <- match.arg(type)

  if (type == "residue_level") {
    f <- system.file("extdata", "Freq_CYS_mutants_res_level.tsv",
                     package = "Cysinger")
    if (f == "") stop("Example data file not found. Is Cysinger installed?")
    # Skip the header line (has fewer tab-separated labels than data columns)
    df <- read.delim(f, header = FALSE, sep = "\t", stringsAsFactors = FALSE,
                     skip = 1)
    # Determine number of bin columns
    n_bins <- ncol(df) - 3
    names(df) <- c("position", "wt_residue", "mut_residue",
                   paste0("bin_", seq_len(n_bins)))
    return(df)
  }

  if (type == "wt_sequence") {
    f <- system.file("extdata", "ccdB_wt_sequence.fasta",
                     package = "Cysinger")
    if (f == "") stop("WT sequence file not found. Is Cysinger installed?")
    lines <- readLines(f, warn = FALSE)
    seq_lines <- lines[!grepl("^>", lines)]
    return(paste0(seq_lines, collapse = ""))
  }
}

#' Generate Example Normalized Data for Demonstrations
#'
#' @description Creates a ready-to-use normalized dataset from the built-in
#' CcdB sample data. Useful for quick demonstrations and testing
#' visualization functions.
#'
#' @return Data frame with columns: position, wt_residue, mut_residue,
#'   bin columns, enrichment columns, and log2 enrichment columns.
#' @export
generate_example_normalized <- function() {
  cys_data <- load_example_data("residue_level")
  normalize_cys_frequency(cys_data, input_bin = 1, sorted_bins = 2:4)
}

#' Run Full Pipeline
#'
#' @description Convenience wrapper to run the full analysis pipeline.
#' Takes FASTQ files and returns cysteine frequency data.
#'
#' @param fwd_fastq Path to forward FASTQ file.
#' @param rev_fastq Path to reverse FASTQ file.
#' @param wt_seq Character. Wild-type nucleotide sequence.
#' @param mid_tags Character vector of MID tags.
#' @param primers Character vector of primers.
#' @param n_bins Integer. Number of bins.
#' @param min_quality Integer. Quality threshold (default 20).
#' @param min_length Integer. Minimum read length (default 75).
#'
#' @return A list with all pipeline outputs.
#' @export
run_pipeline <- function(fwd_fastq, rev_fastq, wt_seq, mid_tags, primers,
                          n_bins = 4L, min_quality = 20L, min_length = 75L) {
  message("Step 1: Filtering forward reads...")
  fwd_reads <- filter_fastq_reads(fwd_fastq, mid_tags, primers,
                                   direction = "forward")
  message("Step 2: Filtering reverse reads...")
  rev_reads <- filter_fastq_reads(rev_fastq, mid_tags, primers,
                                   direction = "reverse")

  message("Step 3: Demultiplexing...")
  fwd_demux <- demultiplex_reads(fwd_reads, mid_tags, primers)
  rev_demux <- demultiplex_reads(rev_reads, mid_tags, primers)

  message("Step 4: Quality filtering...")
  fwd_qc <- quality_filter_reads(fwd_demux, min_quality, min_length)
  rev_qc <- quality_filter_reads(rev_demux, min_quality, min_length)

  message("Step 5: Converting to FASTA...")
  fwd_fasta <- reads_to_fasta(fwd_qc$reads)
  rev_fasta <- reads_to_fasta(rev_qc$reads)

  message("Step 6: Aligning to reference...")
  fwd_aln <- align_to_reference(fwd_fasta$fasta_sequence, wt_seq,
                                 fwd_fasta$read_id)
  rev_aln <- align_to_reference(rev_fasta$fasta_sequence, wt_seq,
                                 rev_fasta$read_id)

  message("Step 7: Combining paired reads...")
  combined <- combine_paired_reads(fwd_aln, rev_aln,
                                    nchar(wt_seq), wt_seq)

  message("Step 8: Building substitution table...")
  sub_tables <- list()
  for (b in seq_len(n_bins)) {
    bin_reads <- combined$combined  # In production, filter by bin
    sub_tables[[b]] <- build_substitution_table(bin_reads, bin_id = b)
  }
  sub_table <- do.call(rbind, sub_tables)

  message("Step 9: Classifying mutations...")
  classified <- classify_mutations(sub_table)

  message("Step 10: Calculating cysteine frequencies...")
  cys_freq <- calculate_cys_frequency(classified$single, wt_seq, n_bins)
  cys_res <- aggregate_cys_residue_level(cys_freq, n_bins)
  cys_norm <- normalize_cys_frequency(cys_res, input_bin = 1,
                                       sorted_bins = 2:n_bins)

  message("Pipeline complete!")

  list(
    fwd_reads = fwd_reads,
    rev_reads = rev_reads,
    fwd_qc_stats = fwd_qc$stats,
    rev_qc_stats = rev_qc$stats,
    combine_stats = combined$stats,
    substitution_table = sub_table,
    classification = classified,
    cys_frequency_codon = cys_freq,
    cys_frequency_residue = cys_res,
    cys_normalized = cys_norm
  )
}
