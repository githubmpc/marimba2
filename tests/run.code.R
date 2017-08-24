## code for testing on different regions with results
## regions to test 76, 153, 192, 351
## simulate under genetic model
set.seed(123)
load_all()
library(tidyverse)
set.seed(123)
stats <- statistics_1000g(region=351)
truth <- simulate_data(stats, N=500, error=0)
##dat <- longFormat(truth$logr, truth$cn)
ggplot(dat, aes(logr, ..count.., fill=copy_number)) +
  geom_density(alpha=0.5) + xlab("LRR")

m <- gmodel(truth$logr, mp=mcmcParams(burnin=0, iter=1000, thin=1))
fit <- gibbs_genetic(m)
gg_chains(fit$chains, truth)
gg_model(fit, bins=100)
gmodel2 <- startAtTrueValues(truth, gmodel)
gg_model(gmodel2, bins=100)
fit2 <- gibbs_genetic(gmodel2, gparams)
gg_chains(fit2$chains, truth)



##trace(gg_model, browser)
gg_model(fit2$gmodel, bins=100)

gibbs.genetic <- gibbs.cnv.wrapper(N=500, K=3, y=y, thin=2,
                                   iter=1000,
                                   xi=100, burnin=0, ncp=30,
                                   model="Genetic")
#t0 <- proc.time()
gibbs.env <- gibbs.cnv.wrapper(N=500, K=3, y=y, thin=2,
                               iter=2000, xi=100,
                               burnin=0, ncp=30,
                               model="Environmental")
#t0 <- proc.time() - t0
gibbs.intermed <- gibbs.cnv.wrapper(N=500, K=3, y=y, thin=2,
                                    iter=2000, xi=100,
                                    burnin=0, ncp=30,
                                    model="Intermediate")
results.genetic <- gibbs.results(truth,gibbs.genetic, gparams)
results.env <- gibbs.results(truth,gibbs.env, gparams)
results.intermed <- gibbs.results(truth,gibbs.intermed, gparams)

# plots
gg_prob.comparison(results.genetic)
gg_prob.comparison(results.env)
gg_prob.comparison(results.intermed)

gg_inten.comparison1(truth, gibbs.genetic, gparams)
gg_inten.comparison1(truth, gibbs.env, gparams)
gg_inten.comparison1(truth, gibbs.intermed, gparams)

gg_inten.comparison2(truth, gibbs.genetic, gparams)
gg_inten.comparison2(truth, gibbs.env, gparams)
gg_inten.comparison2(truth, gibbs.intermed, gparams)

jpeg("logllcomparisons.components.351.jpg")
par(mfrow=c(3,1))
plot(gibbs.genetic$logll, main = "Genetic Model", type = "l",
     #xlim = c(0, 30000), 
     ylab = "Log Likelihood")
plot(gibbs.env$logll, main = "Environmental Model",
     #ylim = c(50, 250), 
     type = "l",
     ylab = "Log Likelihood")
plot(gibbs.intermed$logll, main = "Intermediate Model",
     #ylim = c(50, 250), 
     type = "l",
     ylab = "Log Likelihood")
dev.off()

jpeg("tau.plots.153.jpg")
par(mfrow=c(3,1))
plot(gibbs.genetic$taus, main = "Genetic Model", type = "l",
     #xlim = c(0, 30000), 
     ylab = "Taus")
plot(gibbs.env$taus, main = "Environmental Model",
     #ylim = c(50, 250), 
     type = "l",
     ylab = "Taus")
plot(gibbs.intermed$taus, main = "Intermediate Model",
     #ylim = c(50, 250), 
     type = "l",
     ylab = "Taus")
dev.off()

# extra code
pm.theta <- colMeans(gibbs.genetic$post.theta)
df <- data.frame(y=as.numeric(sim.truth$y),
                 cn=factor(cn.model),
                 truth=factor(as.numeric(sim.truth$cn)),
                 concordant=cn.model==as.numeric(sim.truth$cn))
df[df$cn != df$truth, ]

p <- cbind(as.numeric(gibbs.cnv$post.cn0), as.numeric(gibbs.cnv$post.cn1), as.numeric(gibbs.cnv$post.cn2))
cn.model <- apply(p, 1, which.max) - 1
cn.trio <- matrix(cn.model, params$N, params$K)
colnames(cn.trio) <- c("m", "f", "o")
tmp <- cn.trio != truth$cn
index <- which(rowSums(tmp) > 0)
df[index, ]
truth$cn[index, ]
cn.trio[index, ]
truth$y[index, ]

# manual debug
y <- truth$y
gparams <- geneticParams(N=81, K=3,
                         iter=1000,
                         thin=2,
                         burnin=0, ncp=100, model="Genetic")
gmodel.k3 <- initialize_gmodel(y, gparams, 1)
params <- gparams
gmodel <- gmodel.k3

