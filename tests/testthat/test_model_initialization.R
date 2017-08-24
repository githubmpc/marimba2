context("model initialization")

test_that("geneticParams", {
  gparams <- geneticParams(N=81, K=3)
  expect_is(gparams, "list")
  if(FALSE){
    x <- MCMCpack::rinvgamma(1e6, a, b)
    expect_equal(m, mean(x), tolerance=0.001)
    expect_equal(s, sd(x), tolerance=0.02)
  }
})

test_that("initialize_chains", {
  path <- system.file("extdata", package="marimba")
  truth <- readRDS(file.path(path, "simulation_easy.rds"))
  gparams <- geneticParams(N=81, K=3)
  N <- gparams$N
  comp <- 1
  gmod <- initialize_gmodel(truth$y, gparams, comp)
  K <- length(gmod$theta)
  ch <- initialize_chains(gmod, gparams)
  C <- list()

  for (k in 1:K){
    cn.k <- attributes(table(gmod$cn))$dimnames[[1]][k]
    label <- as.vector(paste0("C",cn.k))
    c <- matrix(0, N, 3)
    colnames(c) <-  c("m", "f", "o")
    C[[length(C)+1]] <- c
    names(C)[k] <- label
  }

  expect_equal(dim(C[[1]]), c(gparams$N, gparams$K))
  expect_equal(length(C), length(gmod$theta))
  expect_is(ch, "list")
})


test_that("initialize_gmodel", {
  path <- system.file("extdata", package="marimba")
  truth <- readRDS(file.path(path, "simulation_easy.rds"))
  set.seed(668899)
  gparams <- geneticParams(N=81, K=3)
  comp <- 1
  ##
  ## gmodel contains parameters for current iteration of the chain
  ##
  # chains <- initialize_chains(gparams)
  ##trace(initialize_gmodel, browser)
  set.seed(149)
  gmodel <- initialize_gmodel(truth$y, gparams, comp)
  table(gmodel$cn, truth$cn)
  expect_true(sum(gmodel$cn != truth$cn) <= 4)
  expect_equal(gmodel$theta, truth$theta, tolerance=0.03)
  expect_equal(gmodel$sigma, truth$sigma, tolerance=0.03)
  cn.parents <- gmodel$cn[, c("m", "f")]
  tab <- table(cn.parents)
  phat <- tab/sum(tab)
  ##
  ## currently, there is no prior on p that supports HWE
  ##
  expect_equal(gmodel$pi, truth$p, tolerance=0.15)
})
