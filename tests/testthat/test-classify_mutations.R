library(testthat)
library(Cysinger)

test_that("classify_mutations categorizes correctly", {
  # Mock substitution table
  mock_sub <- data.frame(
    read_id = c("r1", "r2", "r3", "r4", "r5"),
    bin_id = rep(1, 5),
    n_substitutions = c(0, 1, 2, 3, 5),
    stringsAsFactors = FALSE
  )

  result <- classify_mutations(mock_sub)

  expect_type(result, "list")
  expect_true("wt" %in% names(result))
  expect_true("single" %in% names(result))
  expect_true("summary" %in% names(result))

  # Check correct classification
  expect_equal(nrow(result$wt), 1)
  expect_equal(nrow(result$single), 1)
  expect_equal(nrow(result$double), 1)

  # Summary should have 5 classes
  expect_true(nrow(result$summary) >= 4)
})
