context("Simulating data")

test_that("simulate_data", {
  set.seed(668899)
  ##mendelian.probs <- mendelianProb(epsilon=0)
  mendelian.probs <- gMendelian(tau.one=1, tau.two=0.5, tau.three=0, err=0)
  path <- system.file("extdata", package="marimba")
  vcf.del.df.4 <- read.table(file.path(path, "parameter.p.1000gp.v2.txt"),
                             header=TRUE, sep="\t")
  median.sd.df.4 <- read.table(file.path(path, "medians.sd.txt"),
                               header=TRUE, sep="\t")
  i <- 108 ## region 108 
  p.1000gp <- vcf.del.df.4[i, c(9:11)]
  p.1000gp <- as.numeric(p.1000gp)
  theta.1000gp <- median.sd.df.4[i, c(1:3)]
  theta.1000gp <- as.numeric(theta.1000gp)
  sigma.1000gp <- median.sd.df.4[i, c(4:6)]
  sigma.1000gp <- as.numeric(sigma.1000gp)
  ## custom input as required (optional)
  ##sigma.1000gp[1] <- .08
  ##p.1000gp <- c(.1, .2, .7)
  ##dat <- simulate_data(81, p = p.1000gp,
  ##                     theta = theta.1000gp,
  ##                     sigma = sigma.1000gp,
  ##                     mendelian.probs)
  params <- statistics_1000g(region=108)
  expect_equal(sum(params[, "p"]), 1)
  #params[, "p"] <- c(0.1, 0.2, 0.7)
  #params[1, "sigma"] <- 0.08
  set.seed(668899)
  dat2 <- simulate_data2(params, N=81, error=0)
  ##y1 <- dat$response
  ##y2 <- dat2$y
  ##colnames(y1) <- colnames(y2) <- NULL
  ##expect_identical(y1, y2)
  if(FALSE) {
    saveRDS(dat2, "../../inst/extdata/simulation_easy.rds")
    gg_cnp(dat2)
  }
})
