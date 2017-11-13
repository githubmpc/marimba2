options(error = recover)

library(reshape2)
# library(marimba)
library(gtools)
library(testthat)
library(Hmisc)
library(matrixStats)
library(devtools)
library(coda)
library(tidyverse)
library(magrittr)
library(HWEBayes)
library(coda)
library(mclust)
# load_all("marimba")
setwd("~/Desktop/Chakravarti_Lab/git")
load_all("marimba_two")


p <- c(0.25, 0.5, 0.25)
theta <- c(-4,-1, 2)
sigma <- c(0.3, 0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))

dat2 <- simulate_data(params, N=500, error=0)
y <- dat2$data$log_ratio

gg_cnp(dat2)

gmodel.test<-gmodel2(dat2$data)
gmodel.test<-gmodel(dat2$data)

gibbs.genetic <- gibbs.cnv.call(K=3,
                                   states=0:2,
                                   tau=c(1, 0.5, 0),
                                   xi=c(1.5, 1, 1), 
                                   mu=c(-3, -0.5, 0), ## theoretical
                                   nu=1,
                                   sigma2.0=0.001,
                                   a=c(1, 1, 1),
                                   eta=c(0.5, 0.5),
                                   error=1e-4,
                                   ncp=30,
                                   model="Genetic", dat2)

test.start <- multipleStarts(dat=dat2$data,nstarts=3,iter=10)
test.mcmc <- mcmcList(test.start)

gp=geneticParams()
gp=geneticParams(K=5, states=0:4, xi=c(1.5, 1, 1, 1.5, 1.5), mu=c(-3, -0.5, 0.5, 1.5, 2.5))
gp=geneticParams(K=4, states=0:3, xi=c(1.5, 1, 1, 1.5), mu=c(-3, -0.5, 0.5, 1.5))
gp=geneticParams(K=4, states=1:4, xi=c(1, 1, 1.5, 1.5), mu=c(-0.5, 0.5, 1.5, 2.5))
gp=geneticParams(K=3, states=0:2, xi=c(1.5, 1, 1), mu=c(-3, -0.5, 1))
gp=geneticParams(K=3, states=1:3, xi=c(1.5, 1, 1), mu=c(-0.5, 0.5, 1.5))
gp=geneticParams(K=3, states=2:4, xi=c(1, 1.5, 1.5), mu=c(0.5, 1.5, 2.5))
gp=geneticParams(K=2, states=0:1, xi=c(1.5, 1 ), mu=c(-3, -0.5))
gp=geneticParams(K=2, states=1:2, xi=c(1, 1), mu=c(-0.5, 0.5))
gp=geneticParams(K=2, states=2:3, xi=c(1, 1 ), mu=c(0.5, 1.5))
gp=geneticParams(K=2, states=3:4, xi=c(1.5, 1.5 ), mu=c(1.5, 2.5))

# generate appropriate matrix - must do step
mprob.matrix(tau=c(0.5, 0.5, 0.5), gp=gp)

mp=mcmcParams(burnin=1000, iter=100, thin=1, nstarts=20, max_burnin=3000)
mp=mcmcParams(burnin=20, iter=10, thin=1, nstarts=3, max_burnin=21)
mp=mcmcParams(burnin=20, iter=10, thin=1, nstarts=50, max_burnin=21)

# this set of parameters takes 2.2 hours
mp=mcmcParams(burnin=1000, iter=1000, thin=5, nstarts=5, max_burnin=3000)


mp=mcmcParams(burnin=1500, iter=1000, thin=10, nstarts=10, max_burnin=3500)



start.time <- Sys.time()
gibbs.test <- gibbs(mp, gp, dat2$data)
end.time <- Sys.time()
time.taken <- end.time - start.time

# results
gibbs.select <- selectModels(gibbs.test)
gibbs.test2 <- lapply(gibbs.test, "[", gibbs.select)
gibbs.test2 <- gibbs.test[gibbs.select]
gibbs.diag <- diagnostics(gibbs.test2)
gibbs.unlist <- unlistModels(gibbs.test2)

# results output
current_summary(gibbs.unlist)
gibbs.sum <- posterior_summary(gibbs.unlist)
gibbs.dic <- compute_dic(gibbs.unlist)
gibbs.sum
gg_truth(dat2)
gg_model(gibbs.unlist)
gg_chains(gibbs.unlist, dat2)
posterior_difference(gibbs.sum, dat2)

# multiple model selection mode
data <- dat2$data
gibbs.model <- multiple_models(data, mp)

# when debugging
data <- dat2$data




