library(testthat)
library(Cysinger)

test_that("get_codon_table returns complete table", {
  ct <- get_codon_table()
  expect_length(ct, 64)
  expect_equal(ct[["ATG"]], "M")
  expect_equal(ct[["TGT"]], "C")
  expect_equal(ct[["TGC"]], "C")
  expect_equal(ct[["TAA"]], "*")
  expect_equal(ct[["TAG"]], "*")
  expect_equal(ct[["TGA"]], "*")
})

test_that("translate_sequence works correctly", {
  expect_equal(translate_sequence("ATGTGT"), "MC")
  expect_equal(translate_sequence("ATGTAA"), "M*")
  expect_equal(translate_sequence("ATG"), "M")
  # Incomplete codon should be empty or omitted
  expect_equal(nchar(translate_sequence("AT")), 0)
})

test_that("ccdB_wt_nucleotide returns correct length", {
  wt_nt <- ccdB_wt_nucleotide()
  expect_type(wt_nt, "character")
  expect_equal(nchar(wt_nt), 306)
  # Starts with ATG
 expect_true(startsWith(wt_nt, "ATG"))
})

test_that("ccdB_wt_protein returns correct length", {
  wt_prot <- ccdB_wt_protein()
  expect_type(wt_prot, "character")
  # 101 residues + stop codon
  expect_equal(nchar(wt_prot), 102)
  expect_true(startsWith(wt_prot, "M"))
  expect_true(endsWith(wt_prot, "*"))
})

test_that("reverse_complement_seq is correct", {
  expect_equal(reverse_complement_seq("ATGC"), "GCAT")
  expect_equal(reverse_complement_seq("AAAA"), "TTTT")
  expect_equal(reverse_complement_seq("CCCC"), "GGGG")
  expect_equal(reverse_complement_seq("A"), "T")
  # Should be case-insensitive
  expect_equal(reverse_complement_seq("atgc"), "GCAT")
})
