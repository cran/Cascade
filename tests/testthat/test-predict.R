skip_on_cran()
skip_on_ci()

test_that("predict() returns a 'micropredict' object with coherent slots", {
  skip_if_not_installed("igraph")
  skip_if_not_installed("magic")
  data(M, package = "Cascade")
  data(Net, package = "Cascade")
  # Use a small target set
  targets <- c(5L, 10L)
  pred <- predict(M, Omega = Net, nv = 0, targets = targets)
  expect_s4_class(pred, "micropredict")
  # Slots and classes
  expect_s4_class(pred@network, "network")
  expect_s4_class(pred@microarray_unchanged, "micro_array")
  expect_s4_class(pred@microarray_changed, "micro_array")
  expect_s4_class(pred@microarray_predict, "micro_array")
  expect_equal(pred@nv, 0)
  expect_equal(as.integer(pred@targets), targets)
})
