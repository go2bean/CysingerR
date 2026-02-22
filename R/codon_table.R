#' Codon Translation Table
#'
#' @description Returns the standard genetic code as a named character vector.
#'
#' @return A named character vector mapping codons to amino acids.
#' @export
#' @examples
#' ct <- get_codon_table()
#' ct["ATG"]  # "M"
#' ct["TGT"]  # "C"
get_codon_table <- function() {
  codons <- c(
    "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
    "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "TAT" = "Y", "TAC" = "Y", "TAA" = "*", "TAG" = "*",
    "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
    "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
    "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
    "TGT" = "C", "TGC" = "C", "TGA" = "*", "TGG" = "W",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
  )
  return(codons)
}

#' Translate Nucleotide Sequence to Protein
#'
#' @description Translates a nucleotide sequence to amino acid sequence using
#' the standard genetic code.
#'
#' @param nt_seq Character string of the nucleotide sequence (length must be
#'   divisible by 3).
#' @return Character string of the amino acid sequence.
#' @export
#' @examples
#' translate_sequence("ATGCAGTTTAAG")  # "MQFK"
translate_sequence <- function(nt_seq) {
  nt_seq <- toupper(gsub("[^A-Za-z]", "", nt_seq))
  if (nchar(nt_seq) %% 3 != 0) {
    warning("Sequence length is not divisible by 3; truncating to nearest codon.")
    nt_seq <- substr(nt_seq, 1, (nchar(nt_seq) %/% 3) * 3)
  }
  ct <- get_codon_table()
  codons <- substring(nt_seq, seq(1, nchar(nt_seq), by = 3),
                      seq(3, nchar(nt_seq), by = 3))
  aas <- ct[codons]
  aas[is.na(aas)] <- "X"
  paste0(aas, collapse = "")
}

#' CcdB Wild-Type Nucleotide Sequence
#'
#' @description Returns the CcdB wild-type nucleotide sequence used as the
#' reference in the original cysteine scanning library.
#'
#' @return Character string (306 nt).
#' @export
ccdB_wt_nucleotide <- function() {
  paste0(
    "ATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAG",
    "AGTGATATTATTGACACGCCGGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGTCTGCTG",
    "TCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATG",
    "ATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGC",
    "CACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGAATATAA"
  )
}

#' CcdB Wild-Type Protein Sequence
#'
#' @description Returns the translated CcdB wild-type protein sequence.
#'
#' @return Character string (101 aa + stop).
#' @export
ccdB_wt_protein <- function() {
  translate_sequence(ccdB_wt_nucleotide())
}
