context("Gibbs updates")

test_that("p_offspring", {
  path <- system.file("extdata", package="marimba")
  truth <- readRDS(file.path(path, "simulation_easy.rds"))
  gparams <- geneticParams(N=81, K=3)
  gmod <- initialize_gmodel(truth$y, gparams)
  theta <- gmod$theta
  m.probs <- gmod$m.probs
  cn.mf <- gmod$cn[, c("m", "f")]

  ## heterozygous father, homozygous mother
  expected <- c(0.5, 0.5, 0)
  p <- p_offspring(matrix(c(1, 0), 1, 2), m.probs, theta)
  expect_equal(sum(p), 1)
  p<-as.numeric(p)
  expect_equal(p, expected, tolerance=0.01)

  ## homozygous father, heterozygous mother
  p <- p_offspring(matrix(c(0, 1), 1, 2), m.probs, theta)
  p<-as.numeric(p)
  expect_equal(p, expected, tolerance=0.01)

  ## homozygous father, homozygous mother
  expected <- c(1, 0, 0)
  p <- p_offspring(matrix(c(0, 0), 1, 2), m.probs, theta)
  p<-as.numeric(p)
  expect_equal(p, expected, tolerance=0.01)

  ## diploid father, diploid mother
  expected <- c(0, 0, 1)
  p <- p_offspring(matrix(c(2, 2), 1, 2), m.probs, theta)
  p<-as.numeric(p)
  expect_equal(p, expected, tolerance=0.01)

  ## diploid father, heterozygous mother
  expected <- c(0, 0.5, 0.5)
  p <- p_offspring(matrix(c(2, 1), 1, 2), m.probs, theta)
  p<-as.numeric(p)
  expect_equal(p, expected, tolerance=0.01)

  ## combine 2 rows from above
  expected <- rbind(c(0, 0.5, 0.5), c(0, 0, 1))
  p <- p_offspring(matrix(c(2, 1, 2, 2), 2, 2, byrow=TRUE), m.probs, theta)
  attributes(expected)<-attributes(p)
  expect_equal(p, expected, tolerance=0.01)
})


test_that("cnProb", {
  path <- system.file("extdata", package="marimba")
  truth <- readRDS(file.path(path, "simulation_easy.rds"))
  gparams <- geneticParams(N=81, K=3)
  set.seed(204)
  comp <- 1
  gmod <- initialize_gmodel(truth$y, gparams, comp)
  theta <- gmod$theta

  ## update offspring
  m.probs <- gmod$m.probs
  cn.mf <- gmod$cn[, c("m", "f")]
  p.o <- p_offspring(cn.mf, m.probs, theta)
  y.o <- gmod$y[, "o"]
  theta <- gmod$theta
  sigma <- gmod$sigma
  ##trace(cnProb, browser)
  cn.prob.o <- cnProb(p.o, y.o, theta, sigma)
  cn.prob.o <- head(cn.prob.o)
  set.seed(498)
  z.offspring <- Hmisc::rMultinom(cn.prob.o, 1)[, 1] - 1
  expect_equivalent(z.offspring, c(0, 1, 2, 2, 1, 2))

  ##update parents
  y.m <- gmod$y[, "m"]
  y.f <- gmod$y[, "f"]
  y <- c(y.m, y.f)
  p <- gmod$p
  cn.prob <- cnProb2(p, y.m, theta, sigma)
  cn.prob2 <- cnProb2(p, y, theta, sigma)
  expect_true(all.equal(cn.prob2[1:81, ], cn.prob))
  cn.f <- head(Hmisc::rMultinom(cn.prob2, 1)[, 1] -1, 3)
  expect_equal(c(1, 2, 2), cn.f)
})

test_that("update_cn", {
  ##path <- system.file("extdata", package="marimba")
  ##truth <- readRDS(file.path(path, "simulation_easy.rds"))
  params <- statistics_1000g(region=192)
  truth <- simulate_data(params, N=81, error=0)
  ##
  ##
  ##
  if(FALSE){
    ##gg_model(truth)
  }
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0)
  set.seed(204)
  ##
  ## check model initialization.  Thetas are far enough apart that we should be able to initialize z better
  ##
  ##trace(initialize_gmodel, browser)
  gmod <- initialize_gmodel2(truth$y, gparams)
  N <- gparams$N
  pi <- gmod$pi
  ##gmod2 <- update_mendel(gmod)
  ##
  ## For one iteration, make unit test lenient OR check that the CN estimates
  ## are similar to the starting values
  ##
  table(gmod$cn, truth$cn)
  ##trace(update_cn, browser)
  cn <- update_cn(gmod, gparams)
  table(cn, truth$cn)
  expect_true(sum(cn != gmod$cn) < 20)
  ##
  ## RS: I wouldn't expect this to be true after one iteration and with many
  ## 'mistakes' in the initial values
  ## expect_true(sum(cn != dat2$cn) < 6)
  ## expect_true(sum(cn[, 1:2] != dat2$cn[, 1:2]) < 6)
  ## expect_true(sum(cn[, 3] != dat2$cn[, 3]) < 6)
})

test_that("update_sigma", {
  params <- statistics_1000g(region=192)
  truth <- simulate_data2(params, N=81, error=0)
  gparams <- geneticParams(N=81, K=3,
                          iter=1000,
                          thin=2,
                          burnin=0)
  gparams$nu <- 1
  gparams$sigma2.0 <- 0.001
  set.seed(204)
  comp <- 1
  gmodel <- initialize_gmodel(truth$y, gparams, comp)
  gmodel$sigma <- truth$sigma
  gmodel$cn <- truth$cn
  gmodel$pi <- truth$p
  theta <- gmodel$theta <- truth$theta
  sigma <- update_sigma(params=gparams,
                        ns=n(gmodel),
                        theta=theta,
                        ymns=ybar(gmodel, gparams),
                        yvars=yvar(gmodel, gparams))
  expect_equivalent(sigma, truth$sigma, tolerance=0.05)
})

test_that("update_cn.parents", {
  params <- statistics_1000g(region=192)
  set.seed(123)
  truth <- simulate_data2(params, N=81, error=0)
  if(FALSE){
    y <- truth %$% as.tibble(.$y) %>%
      mutate(trio=paste0("trio", 1:nrow(.)))
    z <- truth %$% as.tibble(.$cn) %>%
      mutate(trio=paste0("trio", 1:nrow(.)))
    y <- y %>%
      gather(key="family_member", value="logR", -trio)
    cn <- z %>%
      gather(key="family_member", value="cn", -trio) %>%
      mutate(cn=factor(cn, levels=0:2))
    dat <- left_join(y, cn)
    ggplot(dat, aes(cn, logR)) +
      geom_jitter(width=0.2) + xlab("latent copy number")

  }
  
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0)
  set.seed(123)
  comp <- 1
  gmod <- initialize_gmodel(truth$y, gparams, comp)
  ## initialize z to true values
  gmod$cn <- truth$cn
  gmod$sigma <- truth$sigma
  gmod$theta <- truth$theta
  gmod$pi <- truth$p

  N <- gparams$N
  pi <- gmod$pi
  gmod <- update_mendel(gmod)
  theta <- gmod$theta
  sigma <- gmod$sigma
  y <- gmod$y
  cn.prob.m <- cnProb2(pi, y[, "m"], theta, sigma)
  expect_true(all.equal(rowSums(cn.prob.m), rep(1, nrow(cn.prob.m))))
  cn.prob.f <- cnProb2(pi, y[, "f"], theta, sigma)
  expect_true(all.equal(rowSums(cn.prob.f), rep(1, nrow(cn.prob.f))))
  table(gmod$cn, truth$cn)

  cn.mf <- update_cn.parents(cn.prob.m, cn.prob.f, gmod, gparams)

  ## repeat a large number times and look at probability of the CN assignments
  z0 <- matrix(0, nrow(cn.mf), 2)
  z1 <- matrix(0, nrow(cn.mf), 2)
  z2 <- matrix(0, nrow(cn.mf), 2)
  for(i in 1:500){
    cn.mf <- update_cn.parents(cn.prob.m, cn.prob.f, gmod, gparams)
    z0 <- z0 + (cn.mf==0 * 1L)
    z1 <- z1 + (cn.mf==1 * 1L)
    z2 <- z2 + (cn.mf==2 * 1L)
  }
  prob.cn0 <- z0/500
  prob.cn1 <- z1/500
  prob.cn2 <- z2/500


  ## 14 changes
  table(cn.mf, truth$cn[, 1:2])
  if(FALSE){
    colnames(cn.mf) <- c("m", "f")
    y <- truth %$% as.tibble(.$y) %>%
      mutate(trio=paste0("trio", 1:nrow(.)))
    z <- truth %$% as.tibble(.$cn) %>%
      mutate(trio=paste0("trio", 1:nrow(.)))
    y <- y %>%
      gather(key="family_member", value="logR", -trio)
    cn <- z %>%
      gather(key="family_member", value="cn", -trio) %>%
      mutate(cn=factor(cn, levels=0:2))
    updated.cn <- as.tibble(cn.mf) %>%
      mutate(trio=paste0("trio", 1:nrow(cn.mf))) %>%
      gather(key="family_member", value="cn", -trio) %>%
      mutate(cn=factor(cn, levels=0:2))
    dat <- left_join(y, cn)
    dat.parents <- dat %>% filter(family_member %in% c("m", "f"))
    dat.parents <- dat.parents %>%
      mutate(not_equal=.$cn != updated.cn$cn)

    probs1 <- as.tibble(prob.cn1) %>%
      set_colnames(c("m", "f")) %>%
      mutate(trio=paste0("trio", 1:nrow(.))) %>%
      gather(key="family_member", value="prob", -trio)

    probs2 <- as.tibble(prob.cn2) %>%
      set_colnames(c("m", "f")) %>%
      mutate(trio=paste0("trio", 1:nrow(.))) %>%
      gather(key="family_member", value="prob", -trio)

    dat.cn1 <- left_join(dat.parents, probs1)
    dat.cn2 <- left_join(dat.parents, probs2)

    ggplot(dat.cn2, aes(cn, logR)) +
      geom_jitter(width=0.2, aes(color=prob)) +
      xlab("latent copy number")
  }

  z0 <- matrix(0, nrow(cn.mf), 3)
  z1 <- matrix(0, nrow(cn.mf), 3)
  z2 <- matrix(0, nrow(cn.mf), 3)
  mendel.prob <- gmod$m.probs
  for(i in 1:500){
    cn.mf <- update_cn.parents(cn.prob.m, cn.prob.f, gmod, gparams)
    cn.o <- update_cn.child(gmod, gparams, cn.mf, mendel.prob)
    cn <- as.tibble(cbind(cn.mf, cn.o))
    z0 <- z0 + (cn==0 * 1L)
    z1 <- z1 + (cn==1 * 1L)
    z2 <- z2 + (cn==2 * 1L)
  }
  prob.cn0 <- z0/500
  prob.cn1 <- z1/500
  prob.cn2 <- z2/500

  if(FALSE){
    y <- truth %$% as.tibble(.$y) %>%
      mutate(trio=paste0("trio", 1:nrow(.)))
    z <- truth %$% as.tibble(.$cn) %>%
      mutate(trio=paste0("trio", 1:nrow(.)))
    y <- y %>%
      gather(key="family_member", value="logR", -trio)
    cn <- z %>%
      gather(key="family_member", value="cn", -trio) %>%
      mutate(cn=factor(cn, levels=0:2))
    dat <- left_join(y, cn)

    probs0 <- as.tibble(prob.cn0) %>%
      set_colnames(c("m", "f", "o")) %>%
      mutate(trio=paste0("trio", 1:nrow(.))) %>%
      gather(key="family_member", value="prob", -trio)

    probs1 <- as.tibble(prob.cn1) %>%
      set_colnames(c("m", "f", "o")) %>%
      mutate(trio=paste0("trio", 1:nrow(.))) %>%
      gather(key="family_member", value="prob", -trio)

    probs2 <- as.tibble(prob.cn2) %>%
      set_colnames(c("m", "f", "o")) %>%
      mutate(trio=paste0("trio", 1:nrow(.))) %>%
      gather(key="family_member", value="prob", -trio)

    dat.cn0 <- left_join(dat, probs0)
    dat.cn1 <- left_join(dat, probs1)
    dat.cn2 <- left_join(dat, probs2)

    ggplot(dat.cn2, aes(cn, logR)) +
      geom_jitter(width=0.2, aes(color=prob)) +
      xlab("latent copy number")

    dat.cn2 %>%
      filter(family_member=="o") %>%
      ggplot(., aes(cn, logR)) +
      geom_jitter(width=0.2, aes(color=prob)) +
      xlab("latent copy number")

    dat.cn1 %>%
      filter(family_member=="o") %>%
      ggplot(., aes(cn, logR)) +
      geom_jitter(width=0.2, aes(color=prob)) +
      xlab("latent copy number")

    dat.cn0 %>%
      filter(family_member=="o") %>%
      ggplot(., aes(cn, logR)) +
      geom_jitter(width=0.2, aes(color=prob)) +
      xlab("latent copy number")
  }


  z0 <- matrix(0, nrow(cn.mf), 3)
  z1 <- matrix(0, nrow(cn.mf), 3)
  z2 <- matrix(0, nrow(cn.mf), 3)
  mendel.prob <- gmod$m.probs
  for(i in 1:500){
    cn <- update_cn(gmod, gparams)
    z0 <- z0 + (cn==0 * 1L)
    z1 <- z1 + (cn==1 * 1L)
    z2 <- z2 + (cn==2 * 1L)
  }
  prob.cn0 <- z0/500
  prob.cn1 <- z1/500
  prob.cn2 <- z2/500
  table(cn, truth$cn)
})

test_that("update_cn.child", {
  params <- statistics_1000g(region=192)
  dat2 <- simulate_data2(params, N=81, error=0)
  truth <- dat2
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0)
  comp <- 1
  gmod <- initialize_gmodel(dat2$y, gparams, comp<-1)
  N <- gparams$N
  gmod <- update_mendel(gmod)
  pi <- gmod$pi
  theta <- gmod$theta
  sigma <- gmod$sigma
  y <- gmod$y
  cn.prob.m <- cnProb2(pi, y[, "m"], theta, sigma)
  cn.prob.f <- cnProb2(pi, y[, "f"], theta, sigma)

  table(gmod$cn, truth$cn)
  cn.mf <- update_cn.parents(cn.prob.m, cn.prob.f, gmod, gparams)

  mendel.prob <- gmod$m.probs

  cn.o <- update_cn.child(gmod, gparams, cn.mf, mendel.prob)
  table(cn.o, truth$cn[,3])
  expect_true(sum(cn.o != truth$cn[,3]) < 10)
})

test_that("balance.cn.all", {
  params <- statistics_1000g(region=192)
  dat2 <- simulate_data2(params, N=81, error=0)
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0)
  comp <- 1
  gmodel <- initialize_gmodel(dat2$y, gparams, comp<-1)
  result <- balance.cn.all(gmodel, gparams)

  cn <- rbinom(100, 1, 0.8)
  gmodel$cn <- cn
  set.seed(123)
  result <- balance.cn.all(gmodel, gparams)
  expect_true(length(missingCnState(result)) == 0)

  cn <- rbinom(100, 1, 0.9999)
  gmodel$cn <- cn
  result <- balance.cn.all(gmodel, gparams)
  expect_true(length(missingCnState(result)) == 0)

  cn <- rbinom(100, 1, 0.000001)
  gmodel$cn <- cn
  result <- balance.cn.all(gmodel, gparams)
  expect_true(length(missingCnState(result)) == 0)

  cn <- rbinom(100, 1, 0.5)
  gmodel$cn <- cn
  result <- balance.cn.all(gmodel, gparams)
  expect_true(length(missingCnState(result)) == 0)
})

test_that("updateTransmissionProb", {
  params <- statistics_1000g(region=192)
  truth <- simulate_data2(params, N=81, error=0)
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0, model="Genetic")
  N <- gparams$N
  K <- gparams$K
  comp <- 1
  gmod <- initialize_gmodel(truth$y, gparams, comp)

  # genetic
  m.probs.correct <- gMendelian(tau.one<-gparams$tau.one, tau.two<-gparams$tau.two, tau.three<-gparams$tau.three, err<-gparams$error)
  gmodel <- updateTransmissionProb(gmod, gparams)
  expect_equal(gmodel$m.probs, m.probs.correct)

  # environmental
  m.probs.correct <- gMendelian(tau.one<-gmod$tau, tau.two<-gmod$tau, tau.three<-gmod$tau, err<-gparams$error)
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0, model="Environmental")
  gmodel <- updateTransmissionProb(gmod, gparams)
  expect_false(isTRUE(all.equal(gmodel$m.probs[,1,1], m.probs.correct[,1,1])))
  expect_false(isTRUE(all.equal(gmodel$m.probs[,2,1], m.probs.correct[,2,1])))
  expect_false(isTRUE(all.equal(gmodel$m.probs[,3,1], m.probs.correct[,3,1])))
  expect_false(isTRUE(all.equal(gmodel$m.probs[,1,2], m.probs.correct[,1,2])))
  expect_false(isTRUE(all.equal(gmodel$m.probs[,2,2], m.probs.correct[,2,2])))
  expect_false(isTRUE(all.equal(gmodel$m.probs[,3,2], m.probs.correct[,3,2])))
  expect_false(isTRUE(all.equal(gmodel$m.probs[,1,3], m.probs.correct[,1,3])))
  expect_false(isTRUE(all.equal(gmodel$m.probs[,2,3], m.probs.correct[,2,3])))
  expect_false(isTRUE(all.equal(gmodel$m.probs[,3,3], m.probs.correct[,3,3])))
  # intermediate
  m.probs.correct <- gMendelian(tau.one<-gparams$tau.one, tau.two<-gmod$tau, tau.three<-gparams$tau.three, err<-gparams$error)
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0, model="Intermediate")
  gmodel <- updateTransmissionProb(gmod, gparams)
  expect_equal(gmodel$m.probs[,1,1], m.probs.correct[,1,1])
  expect_false(isTRUE(all.equal(gmodel$m.probs[,2,1], m.probs.correct[,2,1])))
  expect_equal(gmodel$m.probs[,3,1], m.probs.correct[,3,1])
  expect_false(isTRUE(all.equal(gmodel$m.probs[,1,2], m.probs.correct[,1,2])))
  expect_false(isTRUE(all.equal(gmodel$m.probs[,2,2], m.probs.correct[,2,2])))
  expect_false(isTRUE(all.equal(gmodel$m.probs[,3,2], m.probs.correct[,3,2])))
  expect_equal(gmodel$m.probs[,1,3], m.probs.correct[,1,3])
  expect_false(isTRUE(all.equal(gmodel$m.probs[,2,2], m.probs.correct[,2,3])))
  expect_equal(gmodel$m.probs[,3,3], m.probs.correct[,3,3])
})

# will probably tidy and merge into one update tau function function
test_that("update_tau.env", {
  set.seed(123)
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

  for (i in 1:500) {
    tau <- update_tau.env(gmod, gparams)
  }

  # calculate true taus
  parents.count <- n_parents(truth)
  child.count <- n_child(truth)

  numer.env <- child.count[[2]] + 2 * child.count[[3]]
  denom.env <- parents.count[[2]] + 2 * parents.count[[3]]
  tau.env <- numer.env / denom.env

  expect_equal(tau, tau.env, tolerance=0.15)
})

test_that("update_tau.intermed", {
  set.seed(123)
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

  for (i in 1:500) {
    tau <- update_tau.env(gmod, gparams)
  }

  # calculate true taus
  parents.count <- n_parents(truth)
  child.count <- n_child(truth)

  numer.intermed <- child.count[[2]]
  denom.intermed <- parents.count[[2]]
  tau.intermed <- numer.intermed / denom.intermed

  expect_equal(tau, tau.intermed, tolerance=0.15)
})

test_that("update_pi", {
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
  p <- update_pi(gmod, gparams)
  expect_equivalent(p, gmod$pi, tolerance=0.15)
})

test_that("update_mendel", {
  set.seed(123)
  params <- statistics_1000g(region=192)
  dat2 <- simulate_data2(params, N=81, error=0)
  gparams <- geneticParams(N=81, K=3,
                           iter=1000,
                           thin=2,
                           burnin=0)
  comp <- 1
  gmodel <- initialize_gmodel(dat2$y, gparams, comp<-1)

  gmodel <- update_mendel(gmodel)

  expect_equal(dim(gmodel$m.probs), c(length(gmodel$theta), length(gmodel$theta),length(gmodel$theta)))
})
