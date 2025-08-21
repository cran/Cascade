test_that("micro_array basic structure is consistent", {
  data(M, package = "Cascade")
  expect_s4_class(M, "micro_array")
  expect_true(is.matrix(M@microarray))
  # Columns must equal time points * subjects (validity rule)
  expect_equal(ncol(M@microarray), length(M@time) * M@subject)
  # Basic slot sanity
  expect_true(length(M@time) >= 1)
  expect_true(M@subject >= 1)
})

test_that("printing and showing micro_array work", {
  data(M, package = "Cascade")
#  expect_no_error(show(M))
  expect_no_error(print(M))
})
