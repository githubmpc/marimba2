context("Simple data summaries")

test_that("ybar", {
  params <- statistics_1000g(region=192)
  dat2 <- simulate_data2(params, N=81, error=0)
  truth <- dat2
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0)
  comp <- 1
  gmod <- initialize_gmodel(truth$y, gparams, comp)

  ymns <- ybar(gmod, gparams)

  expect_true(is.numeric(ymns))
  expect_true(!any(is.na(ymns)))

  gmod$y[10,1] <- NA
  ymns <- ybar(gmod, gparams)
  expect_true(is.numeric(ymns))
  expect_true(!any(is.na(ymns)))

  gmod$cn[20,1] <- NA
  ymns <- ybar(gmod, gparams)
  expect_true(is.numeric(ymns))
  expect_true(!any(is.na(ymns)))
})

test_that("yvar", {
  params <- statistics_1000g(region=192)
  dat2 <- simulate_data2(params, N=81, error=0)
  truth <- dat2
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0)
  N <- gparams$N
  K <- gparams$K
  comp <- 1
  gmod <- initialize_gmodel(truth$y, gparams, comp)

  yvars <- yvar(gmod, gparams)
  expect_true(is.numeric(yvars))
  expect_true(!any(is.na(yvars)))

  gmod$cn[20,1] <- NA
  yvars <- yvar(gmod, gparams)
  expect_true(is.numeric(yvars))
  expect_true(!any(is.na(yvars)))

  cn.mat <- matrix(rep.int(1,N*K), nrow = N, ncol = K)
  gmod$cn <- cn.mat
  yvars <- yvar(gmod, gparams)
  expect_true(is.numeric(yvars))
  expect_true(!any(is.na(yvars)))
})
