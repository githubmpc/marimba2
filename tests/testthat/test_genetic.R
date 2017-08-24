##for running locally
if(FALSE){
  library(reshape2)
  library(marimba)
  library(gtools)
  library(testthat)
  library(Hmisc)
  library(matrixStats)
  library(devtools)
  library(coda)
  library(tidyverse)
  library(magrittr)
  load_all("marimba")
}

context("Genetic model for inheritance")

## The recommended method for Stan is to run multiple Markov chains, initialized
## randomly with a diffuse set of initial parameter values, discard the
## warmup/adaptation samples, then split the remainder of each chain in half and
## compute the potential scale reduction statistic, R (Gelman and Rubin, 1992).
## If the result is not enough effective samples, double the number of
## iterations and start again, including rerunning warmup and everything.

test_starts("starts", {
  ##
  ## initialize theta, sigma, p from priors
  ## initial z is from the full conditional
  ##
  set.seed(123)
  stats <- tibble(p=c(0.1, 0.4, 0.5),
                  theta=c(-3, -1, 0),
                  sigma=c(0.2, 0.1, 0.1))
  truth <- simulate_data(stats, N=100)
  gp <- geneticParams()
  mp <- mcmcParams(iter=5, burnin=1, thin=1,
                   nstarts=5, max_burnin=10)
  trace(gibbs, browser)
  fit <- gibbs(mp, gp, truth$data)


  mp <- mcmcParams(iter=100, burnin=1000, thin=1, nstarts=50, max_burnin=10)
  model.list <- gibbs(mp, gp, truth$data)
  diagnostics(model.list)

  current <- .init2(gp)
  current$data <- truth$data
  model <- list(gp=gp,
                mp=mcmcParams(burnin=100, iter=100, thin=1),
                current=current)
  current$z <- update_cn(model)$z



  start.list <- replicate(50, .init2(gp), simplify=FALSE)
  chains <- initialize_chains(states=c(0:2),
                              N=nrow(truth$data),
                              K=3,
                              S=mp$iter+1L)
  fit <- vector("list", nstarts)
  for(i in seq_along(start.list)){
    cat(".")
    curr <- start.list[[i]]
    curr$data <- truth$data
    model <- list(gp=gp, mp=mp, current=curr,
                  chains=chains)
    fit[[i]] <- gibbs_genetic(model)
  }
  mlist <- mcmcList(fit)
  neff <- effectiveSize(mlist)
  r <- gelman.diag(mlist[, -9])
  if(r$mpsrf - 1 < 0.05 && all(neff > 500)){
    ## combine chains
  } else {
    start.list <- replicate(50, .init2(gp), simplify=FALSE)
    ## repeat with twice burnin
    fit <- vector("list", 50)
    mp$burnin <- mp$burnin*2
    for(i in seq_along(start.list)){
      cat(".")
      curr <- start.list[[i]]
      curr$data <- truth$data
      model <- list(gp=gp, mp=mp, current=curr,
                    chains=chains)
      tmp <- gibbs_genetic(model)
      fit[[i]] <- label_switch(tmp)
    }
  }
})

test_that("s1", {
  set.seed(123)
  stats <- tibble(p=c(0.1, 0.4, 0.5),
                  theta=c(-3, -1, 0),
                  sigma=c(0.2, 0.1, 0.1))
  truth <- simulate_data(stats, N=100)
  m <- gmodel(truth$data,
              mp=mcmcParams(burnin=0, iter=100, thin=1))
  dat <- truth$data
  m2 <- gibbs_genetic(m)
  posterior_summary(m2)
  pd <- posterior_difference(posterior_summary(m2), truth)
  expect_identical(nrow(pd$cn.diff), 0L)
  expect_true(all(pd$params$diff < 0.025))
})

test_that("s2", {
  library(purrr)
  set.seed(204)
  stats <- statistics_1000g(region=192)
  truth <- simulate_data(stats, N=100, error=0)

  m <- gmodel(truth$data,
              mp=mcmcParams(burnin=0, iter=100, thin=1))

  if(FALSE){
    gg_truth(truth)
  }
  ## gg_truth(truth)
  mlist <- multipleStarts(dat=truth$data, nstarts=50)
  mp <- list(iter=1000, burnin=500, thin=25)
  mlist2 <- map(mlist, setMcmcParams, mp=mp)
  mlist3 <- map(mlist2, gibbs_genetic)
  mlist3 <- readRDS("mlist3.rds")

  plist <- map(mlist3, posterior_summary)
  pdiff <- map(plist, posterior_difference, truth=truth)
  cn <- do.call(cbind, map(plist, function(x) x$data$copy_number))

  i <- which.max(ll)
  m <- mlist[[i]]
  posterior_difference(current_summary(m), truth)
  m$mp$iter <- 1000
  m$chains <- initialize_chains(states=m$gp$states, N=nrow(m$current$data),
                                )
  m2 <- gibbs_genetic(m)

  m2 <- gibbs_genetic(m)
  posterior_summary(m2)
  pd <- posterior_difference(posterior_summary(m2), truth)
  expect_identical(nrow(pd$cn.diff), 0L)
  expect_true(all(pd$params$diff < 0.025))
})

test_that("region192_model", {
##  mother <- mprob[, , 2] %>% as.tibble %>%
##    gather(key="mother_gt", value="p") %>%
##    mutate(mother_gt=gsub("M_", "", mother_gt))
  set.seed(204)
  stats <- statistics_1000g(region=192)
  truth <- simulate_data(stats, N=81, error=0)
  ##truth$y <- truth$y - median(truth$y)
  ## the theta and sigma chains are highly correlated
  m <- gmodel(truth$data,
              mp=mcmcParams(burnin=250, iter=500, thin=5))
  fit <- gibbs_genetic(m)

  ##
  ## Poor mixing
  if(FALSE)
    gg_chains(fit$chains, truth)
  truez <- truth$cn
  z <- fit$current$z %>% select(c("m", "f", "o")) %>%
    unlist(use.names=FALSE)
  truez <- truth$cn %>% select(c("m", "f", "o")) %>%
    unlist(use.names=FALSE)
  ##fit$current$z
  expect_true(sum(z != truez) < 40)
  cn.map <- map_cn2(fit)
  expect_true(sum(cn.map != truth$cn) <= 23)
  ##
  ## Note:  Above suggests we can do better with more iterations and more thinning.  See effective size:
  ##
  coda::effectiveSize(fit$chains$p)
  pm.theta <- colMeans(fit$chains$theta)
  expect_equivalent(pm.theta, truth$theta, tolerance=0.2)
})

test_that("gmodel_oneiter", {
  path <- system.file("extdata", package="marimba")
  truth <- readRDS(file.path(path, "simulation_easy.rds"))
  gparams <- geneticParams(N=81, K=3)
  set.seed(204)
  comp <- 1
  gmod <- initialize_gmodel(truth$y, gparams, comp)
  N <- gparams$N
  pi <- gmod$pi
  y <- gmod$y
  sigma <- gmod$sigma
  theta <- gmod$theta
  K <- length(theta)

  gmodel <- gmodel_oneiter(gmod, gparams)
  expect_true(sum(gmodel$cn != truth$cn) < 10)
  expect_equal(gmodel$theta, truth$theta, tolerance=0.03)
  expect_equal(gmodel$sigma, truth$sigma, tolerance=0.1)
  expect_equal(gmodel$pi, truth$p, tolerance=0.15)
})

test_that("gibbs_genetic1", {
  params <- statistics_1000g(region=192)
  truth <- simulate_data2(params, N=81, error=0)
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0)
  set.seed(194)
  comp <- 1
  gmodel <- initialize_gmodel(truth$y, gparams, comp)
  cn.initial <- gmodel$cn
  ##trace(gibbs_genetic, browser)
  ##options(error=utils::recover)
  fit <- gibbs_genetic(gmodel, gparams)
  cn.last <- fit$gmodel$cn
  expect_true(!identical(cn.initial, cn.last))
  if(FALSE) {
    gg_cnp(fit$gmodel)
    table(cn.last,  truth$cn)
  }
})

test_that("gibbs_genetic2", {
  path <- system.file("extdata", package="marimba")
  truth <- readRDS(file.path(path, "simulation_easy.rds"))
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0)
  gparams$xi <- 1000
  set.seed(194)
  comp <- 1
  gmodel <- initialize_gmodel(truth$y, gparams, comp)
  fit <- gibbs_genetic(gmodel, gparams)
  S <- gparams$iter + 1
  ch <- fit$chains
  if(FALSE)
    gg_chains(ch)

  df <- data.frame(loglik=ch$logll,
                   iter=seq_len(gparams$iter)*gparams$thin)
  if(FALSE){
    library(ggplot2)
    ggplot(df, aes(iter, loglik)) + geom_line(color="steelblue") +
      geom_hline(yintercept=truth$logll, color="black", linetype="dashed")
  }

  ## compare sigmas
  sigma <- colMeans(ch$sigma)
  names(sigma)<-names(truth$sigma)
  ##empirical.sd <- sapply(split(truth$y, truth$cn), sd)
  expect_equal(sigma, truth$sigma, tolerance=0.11)

  ## compare mix probs
  p <- ch$pi
  p.mns <- as.numeric(colMeans(p))
  cn.truth <- truth$cn
  empirical.p <- as.numeric(table(cn.truth[, 1:2])/(81*2))
  expect_equal(p.mns, empirical.p, tolerance=0.02)

  ## compare thetas
  th <- colMeans(ch$theta)
  names(th)<-names(truth$theta)
  expect_equal(th, truth$theta, tolerance=0.02)
})
