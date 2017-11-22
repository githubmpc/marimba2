# Parameter settings
# CN0 -3
# CN1 -1
# CN2 0.1
# CN3 1.5
# CN4 3.5
# var - 2, 1.5, 1, 1.5, 2


# plot
dat2 <- simulate_data(params, N=500, error=0)
y <- dat2$data$log_ratio
gg_cnp2(dat2)

# multi-allele sim
dat2 <- simulate_data_multi(params, N=500, error=0, gp)
# y <- dat2$data$log_ratio
# gg_cnp2(dat2)

# quick start
mp=mcmcParams(burnin=1000, iter=1000, thin=5, nstarts=5, max_burnin=3000)
data <- dat2
gibbs.model <- multiple_models(data, mp)
# params

# original
p <- c(0.25, 0.5, 0.25)
theta <- c(-4,-1, 2)
sigma <- c(0.3, 0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))

# progressive p = 0.3, q= 0.4, r = 0.3
# p2 + 2pq + (q2 + 2pr) + 2qr + r2 = 1.0 
# K= 5, CN=0,1,2,3,4
p <- c(0.09, 0.24, 0.34, 0.24, 0.09)
theta <- c(-3.5,-1.2, 0.3, 1.7, 4)
sigma <- c(0.3, 0.3, 0.3, 0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))
gp=geneticParams(K=5, states=0:4, xi=c(1.5, 1, 1, 1.5, 1.5), 
                 mu=c(-3, -0.5, 0.5, 1.5, 2.5))

# K = 4, CN = 0,1,2,3
p <- c(0.11, 0.26, 0.37, 0.26)
theta <- c(-3.5,-1.2, 0.3, 1.7)
sigma <- c(0.3, 0.3, 0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))
gp=geneticParams(K=4, states=0:3, xi=c(1.5, 1, 1, 1.5), 
                 mu=c(-3, -0.5, 0.5, 1.5))

# K = 4, CN = 1,2,3,4
p <- c(0.26, 0.37, 0.26, 0.11)
theta <- c(-1.2, 0.3, 1.7, 4)
sigma <- c(0.3, 0.3, 0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))
gp=geneticParams(K=4, states=1:4, xi=c(1, 1, 1.5, 1.5), 
                 mu=c(-0.5, 0.5, 1.5, 2.5))

# p = 0.3
# K= 3, CN=0,1,2
p <- c(0.09, 0.42, 0.49)
theta <- c(-3.5,-1.2, 0.3)
sigma <- c(0.3, 0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))
gp=geneticParams(K=3, states=0:2, xi=c(1.5, 1, 1), 
                 mu=c(-3, -0.5, 1))


# K= 3, CN=1,2,3 (use K=5 probs normalised)
p <- c(0.29, 0.41, 0.30)
theta <- c(-1.2,0.3, 1.7)
sigma <- c(0.3, 0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))
gp=geneticParams(K=3, states=1:3, xi=c(1.5, 1, 1), 
                 mu=c(-0.5, 0.5, 1.5))

# K= 3, CN=2,3,4
p <- c(0.49, 0.42, 0.09)
theta <- c(0.3, 1.7, 4)
sigma <- c(0.3, 0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))
gp=geneticParams(K=3, states=2:4, xi=c(1, 1.5, 1.5), 
                 mu=c(0.5, 1.5, 2.5))

# K= 2, CN=0,1
p <- c(0.18, 0.82)
theta <- c(-3.5, -1.2)
sigma <- c(0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))
gp=geneticParams(K=2, states=0:1, xi=c(1.5, 1 ), 
                 mu=c(-3, -0.5))

# K= 2, CN=1,2
p <- c(0.45, 0.55)
theta <- c(-1.2, 0.3)
sigma <- c(0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))
gp=geneticParams(K=2, states=1:2, xi=c(1, 1), mu=c(-0.5, 0.5))

# K= 2, CN=2,3
p <- c(0.55, 0.45)
theta <- c(0.3, 1.7)
sigma <- c(0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))
gp=geneticParams(K=2, states=2:3, xi=c(1, 1 ), mu=c(0.5, 1.5))

# K= 2, CN=3,4
p <- c(0.82, 0.18)
theta <- c(1.7, 4)
sigma <- c(0.3, 0.3)
params <- data.frame(cbind(p, theta, sigma))
gp=geneticParams(K=2, states=3:4, xi=c(1.5, 1.5 ), mu=c(1.5, 2.5))

