# refer to runcode.onetwo and runcode.three if any bugs running this through in future.

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
# load_all("marimba")
setwd("~/Desktop/Chakravarti_Lab/git")
load_all("marimba2")

setwd("~/Desktop/Chakravarti_Lab/git/marimba2")
##############
####slide 3##
#############
p <- c(0.25, 0.5, 0.25)
theta <- c(-4,-1, 1)
sigma <- c(0.3, 0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))

dat2 <- simulate_data(params, N=500, error=0)
y <- dat2$data$log_ratio

response.df <- melt(dat2$data$log_ratio)
cn.df <- melt(dat2$data$copy_number)

#Visualise
plot.df <- data.frame(logr = response.df$value, cn = cn.df$value)
plot.df$cn <- as.factor(plot.df$cn)
ggplot(plot.df, aes(logr, ..count.., fill = cn)) + 
  geom_density(alpha = .5) #+
  #geom_dotplot(data=mistake.spec.df, aes(LRR.true, fill=factor(CN.true)))

plot.df <- data.frame(logr = response.df$value, cn = cn.df$value)
plot.df$cn <- as.factor(plot.df$cn)
ggplot(plot.df, aes(logr, fill = cn)) + 
  geom_histogram(bins=250) + xlab("LRR")

#############
####slides 4 and 5 simulation replicated 1000x
#############

# simulation
# empty dataframe
parents.count.tbl <- matrix(0, ncol=3, nrow=1000)
child.count.tbl <- matrix(0, ncol=3, nrow=1000)
parents.inten.mean.tbl <- matrix(0, ncol=3, nrow=1000)
child.inten.mean.tbl <- matrix(0, ncol=3, nrow=1000)
parents.inten.var.tbl <- matrix(0, ncol=3, nrow=1000)
child.inten.var.tbl <- matrix(0, ncol=3, nrow=1000)

for (i in 1:1000){
  # simulate data
  dat2 <- simulate_data(params, N=500, error=0)
  y <- dat2$data$log_ratio
  
  # parents
  plot.df <- data.frame(logr = dat2$data$log_ratio, trio = dat2$data$family_member, cn = dat2$data$copy_number)
  plot.df.parents <- plot.df[plot.df$trio == "m" | plot.df$trio == "f",]
  plot.df.child <- plot.df[plot.df$trio == "o",]
  
  # tabulated CN values
  tab.cn.parents <- table(plot.df.parents$cn)
  tab.cn.child <- table(plot.df.child$cn)
  
  # Mean and Sd Intensity values
  tab.inten.mean.parents <- aggregate(plot.df.parents$logr, by=list(Category=plot.df.parents$cn), FUN=mean)
  tab.inten.mean.child <- aggregate(plot.df.child$logr, by=list(Category=plot.df.child$cn), FUN=mean)
  tab.inten.sigma.parents <- aggregate(plot.df.parents$logr, by=list(Category=plot.df.parents$cn), FUN=sd)
  tab.inten.sigma.child <- aggregate(plot.df.child$logr, by=list(Category=plot.df.child$cn), FUN=sd)
  
  # input values
  parents.count.tbl[i,] <- as.numeric(tab.cn.parents)
  child.count.tbl[i,] <- as.numeric(tab.cn.child)
  parents.inten.mean.tbl[i,] <- as.numeric(tab.inten.mean.parents$x)
  child.inten.mean.tbl[i,] <- as.numeric(tab.inten.mean.child$x)
  parents.inten.var.tbl[i,] <- as.numeric(tab.inten.sigma.parents$x)
  child.inten.var.tbl[i,] <- as.numeric(tab.inten.sigma.child$x)
}

# for the 1000 replicates
parents.count.df <- data.frame(parents.count.tbl)
names(parents.count.df) <- c("CN0", "CN1", "CN2")
child.count.df <- data.frame(child.count.tbl)
names(child.count.df) <- c("CN0", "CN1", "CN2")
parents.inten.mean.df <- data.frame(parents.inten.mean.tbl)
names(parents.inten.mean.df) <- c("CN0", "CN1", "CN2")
child.inten.mean.df <- data.frame(child.inten.mean.tbl)
names(child.inten.mean.df) <- c("CN0", "CN1", "CN2")
parents.inten.var.df <- data.frame(parents.inten.var.tbl)
names(parents.inten.var.df) <- c("CN0", "CN1", "CN2")
child.inten.var.df <- data.frame(child.inten.var.tbl)
names(child.inten.var.df) <- c("CN0", "CN1", "CN2")

#need to work out where this comes from....missing plot.parents.cn
plot.parents.cn.2 <- plot.parents.cn
plot.parents.cn.2$value <- plot.parents.cn$value/1000
plot.child.cn.2 <- plot.child.cn
plot.child.cn.2$value <- plot.child.cn$value/500

plot.parents.cn.3 <- plot.parents.cn.2
plot.parents.cn.3$intercept <- ifelse(plot.parents.cn.3$variable=="CN0", 0.25, 0.5)
plot.parents.cn.3$intercept <- ifelse(plot.parents.cn.3$variable=="CN2", 0.25, plot.parents.cn.3$intercept)

ggplot(plot.parents.cn.3, aes(value, colour = variable)) +
  geom_histogram(bins=100) + 
  geom_vline(aes(xintercept=plot.parents.cn.3$intercept), linetype="dotted") +
  facet_grid(.~variable) + 
  ggtitle("Mixture Probabilities (p) for parents across copy number states") +
  xlab("Mixture probabilities (p) per iteration of 500 trios")

plot.child.cn.3 <- plot.child.cn.2
plot.child.cn.3$intercept <- ifelse(plot.child.cn.3$variable=="CN0", 0.25, 0.5)
plot.child.cn.3$intercept <- ifelse(plot.child.cn.3$variable=="CN2", 0.25, plot.child.cn.3$intercept)

ggplot(plot.child.cn.3, aes(value, colour = variable)) +
  geom_histogram(bins=100) + facet_grid(.~variable) + 
  geom_vline(aes(xintercept=plot.child.cn.3$intercept), linetype="dotted") +
  ggtitle("Mixture Probabilities (p) for offspring across copy number states") +
  xlab("Mixture probabilities (p) per iteration of 500 trios")

###############
#####slide 6 HWE check
##############
trio.count <- parents.count.df + child.count.df
trio.paf <- (2*trio.count$CN2 + trio.count$CN1) / (2*(trio.count$CN0+trio.count$CN1+trio.count$CN2))
trio.paf.df <- data.frame(trio.paf)
trio.paf.df$hwe <- "P"
ggplot(trio.paf.df, aes(trio.paf)) + geom_histogram(bins=100) +
  xlab("p") + ggtitle("HWE check for 1000 replicates of 500 simulated trios")

#################################
###slide 8 convergence statistics
##################################
# use region 351 to show poor converging plot
# see run.code.R and run the gibbs sampler there then plot as above
# regions to test 76, 153, 192, 351
# simulate
params <- statistics_1000g(region=133)
#params[, "p"] <- c(0.1, 0.2, 0.7)
#params[1, "sigma"] <- 0.08
set.seed(668899)
dat2 <- simulate_data(params, N=500, error=0)
y <- dat2$y

gparams <- geneticParams(N=500, K=3,
                         iter=1000,
                         thin=5,
                         burnin=250, model="Genetic")

#plot of simulated data
response.df <- melt(dat2$y)
cn.df <- melt(dat2$cn)

plot.df <- data.frame(logr = response.df$value, cn = cn.df$value)
plot.df$cn <- as.factor(plot.df$cn)
ggplot(plot.df, aes(logr, fill = cn)) + 
  geom_density(alpha = .5) + 
  xlab("Log R Ratio (array intensity)") + 
  #ylab("number of samples") +
  ggtitle("Mixture Modelling of CNVs")

gibbs.genetic <- gibbs.cnv.wrapper(N=500, K=2, y=y, thin=5, iter=1000, xi=1000, burnin=250, model="Genetic")
results.genetic <- gibbs.results(dat2,gibbs.genetic)

p.df <- melt(gibbs.genetic$post.pi.child)
colnames(p.df)[2] <- "cn"
#p.df.cn2 <- p.df[p.df$cn==2,]
ggplot(p.df, aes(Var1, value)) +
  geom_line(color="gray")   +  
  geom_hline(data=results.genetic$mixture.prob.tb2, mapping=aes(yintercept=empirical.offspring)) +
  facet_wrap(~cn) +
  ylab("Mixture Probability") +
  xlab("MCMC chain")

################
#######slide 9 histogram of estimated mixture probabilities
###############

mixture.df <- data.frame(gibbs.genetic$post.pi)
theta.df <- data.frame(gibbs.genetic$post.thetas)
names(theta.df) <- c("C0", "C1", "C2")
sigma.df <- data.frame(gibbs.genetic$post.sigmas)
names(sigma.df) <- c("C0", "C1", "C2")

mixture.plot <- melt(mixture.df)
theta.plot <- melt(theta.df)
sigma.plot <- melt(sigma.df)

# HPD interval
library(coda)
mixture.ci <- HPDinterval(mcmc(mixture.df), 0.95)
theta.ci <- HPDinterval(mcmc(theta.df), 0.95)
sigma.ci <- HPDinterval(mcmc(sigma.df), 0.95)

# alternative CI by 0.025 and 0.975
mixture.ci <- apply(gibbs.result$post.pi, 2, quantile, probs = c(0.025, 0.975))
theta.ci <- apply(gibbs.result$post.thetas, 2, quantile, probs = c(0.025, 0.975))
sigma.ci <- apply(gibbs.result$post.sigmas, 2, quantile, probs = c(0.025, 0.975))
mixture.ci <- data.frame(t(mixture.ci))
colnames(mixture.ci) <- c("lower", "upper")
theta.ci <- data.frame(t(theta.ci))
colnames(theta.ci) <- c("lower", "upper")
sigma.ci <- data.frame(t(sigma.ci))
colnames(sigma.ci) <- c("lower", "upper")

# mixture.plot.2 <- mixture.plot[mixture.plot$variable=="C0",]
mixture.plot.2 <- mixture.plot
mixture.plot.2$ci.1 <- ifelse(mixture.plot.2$variable=="C0", mixture.ci[1,1], 0)
mixture.plot.2$ci.1 <- ifelse(mixture.plot.2$variable=="C1", mixture.ci[2,1], mixture.plot.2$ci.1)
mixture.plot.2$ci.1 <- ifelse(mixture.plot.2$variable=="C2", mixture.ci[3,1], mixture.plot.2$ci.1)
mixture.plot.2$ci.2 <- ifelse(mixture.plot.2$variable=="C0", mixture.ci[1,2], 0)
mixture.plot.2$ci.2 <- ifelse(mixture.plot.2$variable=="C1", mixture.ci[2,2], mixture.plot.2$ci.2)
mixture.plot.2$ci.2 <- ifelse(mixture.plot.2$variable=="C2", mixture.ci[3,2], mixture.plot.2$ci.2)

ggplot(mixture.plot.2, aes(value)) + 
  geom_histogram(bins=100)  + 
  facet_grid(. ~ variable, scales = "free") +
  xlab("Mixture Probabilities") +
  #xlim(0.2,0.3) +
  geom_vline(aes(xintercept=mixture.plot.2$ci.1), linetype="dotted") +
  geom_vline(aes(xintercept=mixture.plot.2$ci.2), linetype="dotted") +
  #geom_vline(xintercept=c(mixture.ci[2,1],mixture.ci[2,2]), linetype="dotted") +
  #geom_vline(xintercept=c(mixture.ci[3,1],mixture.ci[3,2]), linetype="dotted") +
  ggtitle("Histogram of Estimated Mixture Probabilities")

###########################
####slide 10 Trace plots
############################
# remember to re-run the gibbs.genetics and results.genetic functions 
# to restore back to original dataset

p.df <- melt(gibbs.genetic$post.pi)
colnames(p.df)[2] <- "cn"
#p.df.cn2 <- p.df[p.df$cn==2,]
ggplot(p.df, aes(Var1, value)) +
  geom_line(color="gray")   +  
  geom_hline(data=results.genetic$mixture.prob.tb2, mapping=aes(yintercept=empirical.parents)) +
  facet_wrap(~cn) +
  ylab("Mixture Probability") +
  xlab("MCMC chain (parents)")

p.df <- melt(gibbs.genetic$post.pi.child)
colnames(p.df)[2] <- "cn"
#p.df.cn2 <- p.df[p.df$cn==2,]
ggplot(p.df, aes(Var1, value)) +
  geom_line(color="gray")   +  
  geom_hline(data=results.genetic$mixture.prob.tb2, mapping=aes(yintercept=empirical.offspring)) +
  facet_wrap(~cn) +
  ylab("Mixture Probability") +
  xlab("MCMC chain (offspring)")

######################
### Running the Gibbs sampler over 50-100 replicates###
######################
iter <- 50

# prep tables
theta.tbl <- matrix(0, ncol=3, nrow=iter)
sigma.tbl <- matrix(0, ncol=3, nrow=iter)
mixture.tbl <- matrix(0, ncol=3, nrow=iter)
mixture.tbl.offspring <- matrix(0, ncol=3, nrow=iter)
parent.misclass.tbl <- matrix(0, ncol=3, nrow=iter)
child.misclass.tbl <- matrix(0, ncol=3, nrow=iter)
misclass.perc <- matrix(0, ncol=2, nrow=iter)
prop.mis.tbl <- matrix(0, ncol=nrow(params)*(ncol(params)+1), nrow=iter)
prop.mis.tbl.q <- matrix(0, ncol=nrow(params)*(ncol(params)+1), nrow=iter)
mistake.typed.tbl <- matrix(0, ncol=8, nrow=1)
ess.tbl <- matrix(0, ncol=nrow(params)*(ncol(params)+1), nrow=iter)
mistake.values.tbl <- matrix(0, ncol=3, nrow=1)

t0 <- proc.time()

options(error = recover) 

for (i in 1:iter){
  # simulate data
  dat2 <- simulate_data(params, N=500, error=0)
  #y <- dat2$y
  y <- matrix(data = dat2$data$log_ratio, nrow = 500, ncol = 3, byrow = T)
  y <- y[, c(2,1,3)]
  colnames(y) <- c("m", "f", "o")
  
  # run model
  # t0 <- proc.time()
  gibbs.genetic <- gibbs.cnv.wrapper(K=3,
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
                                     model="Genetic", y=y)
  results.genetic <- gibbs.results(dat2, gibbs.genetic)
  gibbs.genetic <- ess.check(gibbs.genetic)
  ci.mat <- credible.interval(gibbs.genetic, results.genetic, params)
  ci.mat.q <- credible.interval.quantile(gibbs.genetic, results.genetic, params)
  ess.mat <- ess.compile(results.genetic, params)
  # t0 <- proc.time() - t0
  
  # document results
  theta.tbl[i,] <- apply(gibbs.genetic$post.thetas,2,mean)
  sigma.tbl[i,] <- apply(gibbs.genetic$post.sigmas,2,mean)
  mixture.tbl[i,] <- results.genetic$mixture.prob.tb2$estimated.parents
  mixture.tbl.offspring[i,] <- results.genetic$mixture.prob.tb2$estimated.offspring
  prop.mis.tbl[i,] <- ci.mat
  prop.mis.tbl.q[i,] <- ci.mat.q
  ess.tbl[i,] <- ess.mat
  parent.misclass.tbl[i,] <- results.genetic$parent.misclassification.by.cn
  child.misclass.tbl[i,] <- results.genetic$child.misclassification.by.cn
  misclass.perc[i,1] <- results.genetic$parent.misclassification.perc
  misclass.perc[i,2] <- results.genetic$child.misclassification.perc
  mistake.typed.tbl <- smartbind(mistake.typed.tbl, results.genetic$mistake.typed)
  if (nrow(results.genetic$missed.call.values)!=0){
    # line below changes to match iter
    iteration <- c(1:100)
    iter.no <- data.frame(iteration[i])
    names(iter.no) <- "iteration"
    mistake.values.tbl <- smartbind(mistake.values.tbl, results.genetic$missed.call.values, iter.no)
  }
  
  #stopifnot(mistake.typed.tbl$'0'[i+1]>400)
}

t0 <- proc.time() - t0

colnames(prop.mis.tbl) <- c("pi.0", "pi.1", "pi.2", "pi.child.0", "pi.child.1", "pi.child.2", "theta.0", "theta.1", "theta.2", "sigma.0", "sigma.1", "sigma.2")
colnames(prop.mis.tbl.q) <- c("pi.0", "pi.1", "pi.2", "pi.child.0", "pi.child.1", "pi.child.2", "theta.0", "theta.1", "theta.2", "sigma.0", "sigma.1", "sigma.2")
colnames(ess.tbl) <- c("pi.0", "pi.1", "pi.2", "pi.child.0", "pi.child.1", "pi.child.2", "theta.0", "theta.1", "theta.2", "sigma.0", "sigma.1", "sigma.2")

# adjust for each run for colnames depending on output
# 1 = M, 10 = F, 100 = O, MF = 11, FO = 110, MO = 101, MFO = 111
colnames(mistake.typed.tbl) <- c("Zero", "F", "M", "O", "MO", "FO")
mistake.typed.tbl[is.na(mistake.typed.tbl)] <- 0
mistake.typed.tbl <- mistake.typed.tbl[-1,]
rownames(mistake.typed.tbl) <- c(1:iter)

# adjust to calculate mean per replicate
mistake.typed.tbl$total <-  apply(mistake.typed.tbl,1, sum) - mistake.typed.tbl$Zero
mean(mistake.typed.tbl$total)

# adjust mistake.values.tbl
mistake.values.tbl <- mistake.values.tbl[-1,]

# save outputs for now
write.table(theta.tbl,"theta.10.1000.c4.txt", sep="\t", quote=F)
write.table(sigma.tbl,"sigma.10.1000.c4.txt", sep="\t", quote=F)
write.table(mixture.tbl,"mixture.parents.10.1000.c4.txt", sep="\t", quote=F)
write.table(mixture.tbl.offspring, "mixture.offspring.10.1000.c4.txt")
write.table(parent.misclass.tbl,"misclas.par.10.1000.c4.txt", sep="\t", quote=F)
write.table(child.misclass.tbl,"misclas.child.10.1000.c4.txt", sep="\t", quote=F)
write.table(misclass.perc,"misclas.perc.10.1000.c4.txt", sep="\t", quote=F)
write.table(prop.mis.tbl, "misclass.proportion.ci.c4.txt", sep="\t", quote=F)
write.table(prop.mis.tbl.q, "misclass.proportion.ci.quantile.c4.txt", sep="\t", quote=F)
write.table(mistake.typed.tbl, "mistake.typed.tbl.c4.txt", sep="\t", quote=F)
write.table(ess.tbl, "ess.tbl.c4.txt", sep="\t", quote=F)
write.table(mistake.values.tbl, "mistake.values.tbl.c4.txt", sep="\t", quote=F)

##########################
###error follow up plots
#########################3

gg_chains(gibbs.genetic)
gg_inten.comparison1(dat2,gibbs.genetic)
gg_inten.comparison2(dat2,gibbs.genetic)
gg_theta.comparison(gibbs.genetic)
gg_sigma.comparison(gibbs.genetic)
gg_prob.comparison(gibbs.genetic, results.genetic)
gg_prob.child.comparison(gibbs.genetic, results.genetic)

###################
###results input###
###################

# means
apply(theta.tbl,2,mean)
apply(sigma.tbl,2,mean)
apply(mixture.tbl,2,mean)
apply(mixture.tbl.offspring,2, mean)
apply(parent.misclass.tbl,2,mean)
apply(child.misclass.tbl,2,mean)
apply(prop.mis.tbl,2,mean)
apply(prop.mis.tbl.q,2,mean)
apply(misclass.perc,2,mean)
apply(mistake.typed.tbl,2,sum)

# t test
m <- mistake.typed.tbl$M
f <- mistake.typed.tbl$'F'
o <- mistake.typed.tbl$O
t.test(m,o)
t.test(f,o)
t.test(m,f)

# prep for plots
theta.tbl.df <- data.frame(theta.tbl)
names(theta.tbl.df) <- c("CN0", "CN1", "CN2")
sigma.tbl.df <- data.frame(sigma.tbl)
names(sigma.tbl.df) <- c("CN0", "CN1", "CN2")
mixture.tbl.df <- data.frame(mixture.tbl)
names(mixture.tbl.df) <- c("CN0", "CN1", "CN2")
parent.misclass.tbl.df <- data.frame(parent.misclass.tbl)
names(parent.misclass.tbl.df) <- c("CN0", "CN1", "CN2")
child.misclass.tbl.df <- data.frame(child.misclass.tbl)
names(child.misclass.tbl.df) <- c("CN0", "CN1", "CN2")
misclass.perc.df <- data.frame(misclass.perc)
names(misclass.perc.df) <- c("Parents", "Offspring")

theta.tbl.cn <- melt(theta.tbl.df)
sigma.tbl.cn <- melt(sigma.tbl.df)
plot.mixture.tbl <- melt(mixture.tbl.df)
plot.parent.misclass.tbl <- melt(parent.misclass.tbl.df)
plot.child.misclass.tbl <- melt(child.misclass.tbl.df)
plot.misclass.perc <- melt(misclass.perc.df)

# do HWE p value plot
mixture.tbl.df.hwe <- mixture.tbl.df * 1000

#####################
####slide 12 - accuracy of mixture probability estimation
###############
mixture.prop.ci <- apply(prop.mis.tbl.q[,1:3], 2, sum)
mixture.prop.ci <- mixture.prop.ci/50
plot.mixture.prop.ci <- melt(mixture.prop.ci)
plot.mixture.prop.ci$variable <- c("p0", "p1", "p2")
ggplot(plot.mixture.prop.ci, aes(variable, value)) + geom_point() + coord_flip() +
  xlab("Parameters") + ylab("Proportion within Bayesian Credible Interval")

#########
######slide 13 
###############
#theta.plot.2 <- theta.plot[theta.plot$variable=="C2",]
theta.plot.2 <- theta.plot
theta.plot.2$ci.1 <- ifelse(theta.plot.2$variable=="C0", theta.ci[1,1], 0)
theta.plot.2$ci.1 <- ifelse(theta.plot.2$variable=="C1", theta.ci[2,1], theta.plot.2$ci.1)
theta.plot.2$ci.1 <- ifelse(theta.plot.2$variable=="C2", theta.ci[3,1], theta.plot.2$ci.1)
theta.plot.2$ci.2 <- ifelse(theta.plot.2$variable=="C0", theta.ci[1,2], 0)
theta.plot.2$ci.2 <- ifelse(theta.plot.2$variable=="C1", theta.ci[2,2], theta.plot.2$ci.2)
theta.plot.2$ci.2 <- ifelse(theta.plot.2$variable=="C2", theta.ci[3,2], theta.plot.2$ci.2)

ggplot(theta.plot.2, aes(value)) + 
  geom_histogram(bins=100)  + 
  facet_grid(. ~ variable, scales = "free") +
  xlab("Mean Intensity") +
  #xlab("Estimated copy number state") +
  geom_vline(aes(xintercept=theta.plot.2$ci.1), linetype="dotted") +
  geom_vline(aes(xintercept=theta.plot.2$ci.2), linetype="dotted") +
  ggtitle("Histogram of Estimated Mean Intensity")

# slide 14 bivariate plot of thetas and sigmas

theta.sigma.plot <- cbind(theta.plot, sigma.plot)
theta.sigma.plot <- theta.sigma.plot[,-3]
names(theta.sigma.plot) <- c("variable", "theta", "sigma")

params.df <- data.frame(params)

ggplot(theta.sigma.plot, aes(sigma, theta)) + geom_density_2d() +
ylim(-5,2) + xlim(0.02,0.25) +
geom_point(data=params.df, aes(sigma, theta))

# slide 15
theta.sigma.prop.ci <- apply(prop.mis.tbl.q[,4:9], 2, sum)
theta.sigma.prop.ci <- theta.sigma.prop.ci/50
plot.theta.sigma.prop.ci <- melt(theta.sigma.prop.ci)
plot.theta.sigma.prop.ci$variable <- c("theta.cn0", "theta.cn1", "theta.cn2", "sigma.cn0", "sigma.cn1", "sigma.cn2")
ggplot(plot.theta.sigma.prop.ci, aes(variable, value)) + geom_point() + coord_flip() +
  xlab("Parameters") + ylab("Proportion within Bayesian Credible Interval")

# slide 17 Validating copy number probabilities

gparams <- geneticParams(N=500, K=3,
                         iter=1000,
                         thin=5,
                         burnin=250)
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


## 16 changes
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
 
  # slide 18 misclassification in one simulation
  
  #dat2 <- simulate_data2(params, N=500, error=0)
  # note dat2 is carried over from the loop for the last replicate
  y <- dat2$y
  
  response.df <- melt(dat2$y)
  cn.df <- melt(dat2$cn)
  
 # gibbs.genetic <- gibbs.cnv.wrapper(N=500, K=3, y=y, thin=5, iter=1000, xi=100, burnin=250, model="Genetic")
 # results.genetic <- gibbs.results(dat2, gibbs.genetic)
  
  # choose one, not both
  mistake.spec.df <- data.frame(results.genetic$missed.call.values)
  mistake.spec.df <- data.frame(mistake.values.tbl)
  
  #Visualise mistake on truth distribution
  plot.df <- data.frame(logr = response.df$value, cn = cn.df$value)
  plot.df$cn <- as.factor(plot.df$cn)
  ggplot(plot.df, aes(logr, ..count.., fill = cn)) + 
    geom_density(alpha = .5) +
    geom_dotplot(data=mistake.spec.df, aes(LRR.true, fill=factor(CN.estimated)))
    
# slide 19 misclassification across one simulation
  # adjust for each run for colnames depending on output
  # 1 = M, 10 = F, 100 = O, MF = 11, FO = 110, MO = 101, MFO = 111
  colnames(mistake.typed.tbl) <- c("Zero")
  mistake.typed.tbl <- data.frame(mistake.typed.tbl[-1,])
  
  plot.misclass.type.df <- melt(mistake.typed.tbl)
  plot.misclass.type.df.2 <- plot.misclass.type.df[plot.misclass.type.df$variable!="Zero",]
  plot.misclass.type.df.2$value <- plot.misclass.type.df.2$value/500
  ggplot(plot.misclass.type.df.2, aes(variable, value)) + geom_point() + coord_flip() +
    xlab("Error Type") + ylab("Proportion of misclassified samples") + ylim(0,1)
  
# slide 20 misclassification across 100 replicates

  mistakes.count <- data.frame(apply(mistake.typed.tbl,2,sum))
  #colnames(mistakes.count) <- c("Zero")
  
  plot.mistakes.count <- melt(mistakes.count)
  plot.mistakes.count.2 <- plot.mistakes.count[plot.mistakes.count$variable!="Zero",]
  
  ggplot(plot.mistakes.count.2, aes(variable, value)) + geom_point() + coord_flip() +
    xlab("Error Type") + ylab("Count of Misclassified Samples Across 100 replicates")
  
  # slide 21 framework and parameters for different simulations
  
  # create vector of different thetas
  # theta.inputs<- c(-4, -1,1, -2, -0.5, 1, -1.5, -0.5, 0.25)
  # theta.matrix <- matrix(theta.inputs, nrow=3, ncol=3, byrow=T)
  # sigma.inputs<- c(0.1,0.1,0.1, 0.2, 0.2, 0.2, 0.3, 0.3, 0.3)
  # sigma.matrix <- matrix(sigma.inputs, nrow=3, ncol=3, byrow=T)
  # lrr.matrix <- matrix(0, nrow=1500, ncol=9, byrow=T)
  # cn.matrix <- matrix(0, nrow=1500, ncol=9, byrow=T)
  
  theta.inputs<- c(-4, -1,1, -2, -0.75, 0.5, -1.5, -0.5, 0.25)
  theta.matrix <- matrix(theta.inputs, nrow=3, ncol=3, byrow=T)
  sigma.inputs<- c(0.15,0.15,0.15, 0.2, 0.2, 0.2, 0.25, 0.25, 0.25)
  sigma.matrix <- matrix(sigma.inputs, nrow=3, ncol=3, byrow=T)
  lrr.matrix <- matrix(0, nrow=1500, ncol=9, byrow=T)
  cn.matrix <- matrix(0, nrow=1500, ncol=9, byrow=T)
  
  # substitute for differents ps as appropriate
  #A (p = 0.5, MAF = 0.5) : p <- c(0.25, 0.5, 0.25)
  #B (p = 0.85, MAF = 0.15): p <- c(0.0225, 0.255, 0.7225)
  #C (p = 0.25, MAF = 0.75): p <- c(0.5625, 0.375, 0.0625)
  
  for (i in 1:nrow(theta.matrix)){
    for (j in 1:nrow(sigma.matrix)){
      # p <- c(0.25, 0.5, 0.25)
      p <- c(0.0225, 0.255, 0.7225)
      # p <- c(0.5625, 0.375, 0.0625)
      theta <- theta.matrix[i,]
      sigma <- sigma.matrix[j,]
      params <- cbind(p, theta, sigma)
      
      dat2 <- simulate_data2(params, N=500, error=0)
      y <- dat2$y
      response.df <- melt(dat2$y)
      cn.df <- melt(dat2$cn)
      
      k <- (i-1)*3 + j
      lrr.matrix[,k] <- response.df$value
      cn.matrix[,k] <- cn.df$value
    }
  }
  
  lrr.df <- data.frame(lrr.matrix)
  cn.df <- data.frame(cn.matrix)
  names(lrr.df) <- c("Sim1", "Sim2", "Sim3","Sim4","Sim5","Sim6","Sim7","Sim8","Sim9")
  names(cn.df) <- c("Sim1", "Sim2", "Sim3","Sim4","Sim5","Sim6","Sim7","Sim8","Sim9")
  lrr.df.2 <- melt(lrr.df)
  cn.df.2 <- melt(cn.df)
  plot.df <- cbind(lrr.df.2, cn.df.2)
  plot.df <- plot.df[,c(1,2,4)]
  names(plot.df) <- c("Simulation", "LRR", "CN")
  
  plot.df$cn <- as.factor(plot.df$CN)
  ggplot(plot.df, aes(LRR,  fill = CN, colour=CN)) + 
    facet_wrap(~Simulation, scales="free_y") +
    geom_density(alpha = .5) + 
    xlab("Log R Ratio (array intensity)") + 
    #ylab("number of samples") +
    ggtitle("Simulation Roadmap for MAF = 0.15")
  
  
