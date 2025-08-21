skip_on_cran()
skip_on_ci()

.loess_warn_patterns <- c(
  "singularit",
  "span too small",
  "fewer data values",
  "pseudoinverse used",
  "neighborhood radius",
  "reciprocal condition number"
)

test_that("cutoff() returns p-values and smoothed p-values", {
  skip_if_not_installed("VGAM")
  data(Net, package = "Cascade")
  seq_nv <- seq(0, 0.2, length.out = 5)
  
  res <- muffle_warnings_matching(
    cutoff(Net, sequence = seq_nv),
    patterns = .loess_warn_patterns
  )
  
  expect_type(res, "list")
  expect_true(all(c("p.value","p.value.inter","sequence") %in% names(res)))
  expect_equal(length(res$p.value), length(seq_nv))
  expect_equal(length(res$p.value.inter), length(seq_nv))
  expect_true(all(is.finite(res$p.value)))
  expect_true(all(res$p.value >= 0 & res$p.value <= 1))
})

test_that("compare() returns 5 performance metrics between 0 and 1", {
  data(Net, package = "Cascade")
  data(Net_inf, package = "Cascade")
  nv <- 0
  perf <- Cascade::compare(Net, Net_inf, nv)
  expect_type(perf, "double")
  expect_equal(length(perf), 5L)
  expect_true(all(perf >= 0 & perf <= 1))
})
