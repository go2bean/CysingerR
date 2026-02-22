library(testthat)
library(Cysinger)

test_that("normalize_cys_frequency adds enrichment columns", {
  # Create mock data
  mock_data <- data.frame(
    position = 1:5,
    wt_residue = c("M", "Q", "F", "K", "V"),
    mut_residue = rep("C", 5),
    bin_1 = c(100, 200, 50, 150, 80),
    bin_2 = c(120, 180, 60, 140, 90),
    bin_3 = c(200, 100, 30, 250, 40),
    bin_4 = c(500, 50, 10, 400, 20),
    stringsAsFactors = FALSE
  )

  norm <- normalize_cys_frequency(mock_data, input_bin = 1, sorted_bins = 2:4)

  # Should have enrichment columns
  expect_true("enrichment_bin_2" %in% names(norm))
  expect_true("enrichment_bin_3" %in% names(norm))
  expect_true("enrichment_bin_4" %in% names(norm))

  # Should have log2 enrichment columns
  expect_true("log2_enrichment_bin_2" %in% names(norm))
  expect_true("log2_enrichment_bin_4" %in% names(norm))

  # Enrichment should be positive
  expect_true(all(norm$enrichment_bin_2 >= 0))
  expect_true(all(norm$enrichment_bin_4 >= 0))
})

test_that("load_example_data returns correct structure", {
  # This test will only work when package is installed
  skip_if_not_installed("Cysinger")

  cys_data <- load_example_data("residue_level")
  expect_s3_class(cys_data, "data.frame")
  expect_true("position" %in% names(cys_data))
  expect_true("wt_residue" %in% names(cys_data))
  expect_true("bin_1" %in% names(cys_data))
  expect_equal(nrow(cys_data), 101)
})
