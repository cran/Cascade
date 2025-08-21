test_that("network dataset loads and has expected shape", {
  data(network, package = "Cascade")
  expect_s4_class(network, "network")
  expect_true(is.matrix(network@network))
  expect_equal(nrow(network@network), ncol(network@network))
})

test_that("position() returns a 2-column coordinate matrix", {
  skip_if_not_installed("igraph")
  data(network, package = "Cascade")
  pos <- Cascade::position(network, nv = 0)
  expect_true(is.matrix(pos))
  expect_equal(ncol(pos), 3)
  expect_equal(nrow(pos), nrow(network@network) - length(which(rowSums(abs(network@network) > 0) + colSums(abs(network@network) > 0) == 0)))
})

test_that("analyze_network() returns expected columns", {
  skip_if_not_installed("tnet")
  data(network, package = "Cascade")
  ana <- analyze_network(network, nv = 0)
  expect_s3_class(ana, "data.frame")
  expect_true(all(c("node","betweenness","degree","output","closeness") %in% colnames(ana)))
  expect_equal(nrow(ana), nrow(network@network))
})

test_that("geneNeighborhood() returns neighborhoods when graph=FALSE", {
  skip_if_not_installed("igraph")
  data(network, package = "Cascade")
  set.seed(123)
  tgt <- sample(seq_len(nrow(network@network)), size = 1)
  K <- geneNeighborhood(network, targets = tgt, nv = 0, order = 2, graph = FALSE)
  expect_type(K, "list")
  expect_equal(length(K), 2)
})
