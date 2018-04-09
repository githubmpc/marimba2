context("Simulating data")

test_that("simulate_data_multi", {
  set.seed(98765)
  ##mendelian.probs <- mendelianProb(epsilon=0)
  p <- c(0.09, 0.24, 0.34, 0.24, 0.09)
  theta <- c(-3.5,-1.2, 0.3, 1.7, 4)
  sigma <- c(0.2, 0.2, 0.2, 0.2, 0.2)
  params <- data.frame(cbind(p, theta, sigma))
  gp<-geneticParams(K=5, states=0:4, xi=c(0.2, 0.2, 0.2, 0.2, 0.2), 
                    mu=c(-3.5, -1.2, 0.3, 1.7, 4))
  #mendelian.probs <- gMendelian.multi()
  expect_equal(sum(params[, "p"]), 1)
  
  dat2 <- simulate_data_multi(params, N=500, error=0, gp)
  expect_that(length(levels(factor(dat2$data$copy_number))) == length(theta), is_true())
  expect_equal(p, dat2$params$p, tolerance=0.05)
  expect_equal(theta, dat2$params$theta, tolerance=0.03)
  expect_equal(sigma, dat2$params$sigma, tolerance=0.03)
  
  p <- c(0.11, 0.26, 0.37, 0.26)
  theta <- c(-3.5,-1.2, 0.3, 1.7)
  sigma <- c(0.3, 0.3, 0.3, 0.3)
  params <- data.frame(cbind(p, theta, sigma))
  gp=geneticParams(K=4, states=0:3, xi=c(0.3, 0.3, 0.3, 0.3), 
                   mu=c(-3.5, -1.2, 0.3, 1.7))
  
  dat2 <- simulate_data_multi(params, N=500, error=0, gp)
  expect_that(length(levels(factor(dat2$data$copy_number))) == length(theta), is_true())
  expect_equal(p, dat2$params$p, tolerance=0.05)
  expect_equal(theta, dat2$params$theta, tolerance=0.03)
  expect_equal(sigma, dat2$params$sigma, tolerance=0.03)
  
  p <- c(0.49, 0.42, 0.09)
  theta <- c(0.3, 1.7, 4)
  sigma <- c(0.2, 0.2, 0.2)
  params <- data.frame(cbind(p, theta, sigma))
  gp=geneticParams(K=3, states=2:4, xi=c(0.2, 0.2, 0.2), 
                   mu=c(0.3, 1.7, 4))
  
  dat2 <- simulate_data_multi(params, N=500, error=0, gp)
  expect_that(length(levels(factor(dat2$data$copy_number))) == length(theta), is_true())
  expect_equal(p, dat2$params$p, tolerance=0.05)
  expect_equal(theta, dat2$params$theta, tolerance=0.03)
  expect_equal(sigma, dat2$params$sigma, tolerance=0.03)
  path <- system.file("extdata", package="marimba2")
  saveRDS(dat2, paste0(path,"/simulation_test.rds"))
  
  p <- c(0.55, 0.45)
  theta <- c(0.3, 1.7)
  sigma <- c(0.3, 0.3)
  params <- data.frame(cbind(p, theta, sigma))
  gp=geneticParams(K=2, states=2:3, xi=c(0.3, 0.3), mu=c(0.3, 1.7))
  
  dat2 <- simulate_data_multi(params, N=500, error=0, gp)
  expect_that(length(levels(factor(dat2$data$copy_number))) == length(theta), is_true())
  expect_equal(p, dat2$params$p, tolerance=0.05)
  expect_equal(theta, dat2$params$theta, tolerance=0.03)
  expect_equal(sigma, dat2$params$sigma, tolerance=0.03)
  
})
