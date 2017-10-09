# Simulation testing for p=0.5

# A1 - p = 0.5
p <- c(0.25, 0.5, 0.25)
theta <- c(-4,-0.5, 0.5)
sigma <- c(0.2, 0.1, 0.05)
params <- cbind(p, theta, sigma)

# A2 - increase sigmas, little overlap
p <- c(0.25, 0.5, 0.25)
theta <- c(-4,-1, 1)
sigma <- c(0.3, 0.3, 0.3)
params <- cbind(p, theta, sigma)

# A3 - closer thetas, little overlap
p <- c(0.25, 0.5, 0.25)
theta <- c(-2,-0.75, 0.5)
sigma <- c(0.2, 0.2, 0.2)
params <- cbind(p, theta, sigma)

# A4 - moderate overlap
p <- c(0.25, 0.5, 0.25)
theta <- c(-1,-0.3, 0.25)
sigma <- c(0.15, 0.15, 0.15)
params <- cbind(p, theta, sigma)

# B1 - p = 0.85
p <- c(0.0225, 0.255, 0.7225)
theta <- c(-4,-0.5, 0.5)
sigma <- c(0.2, 0.1, 0.05)
params <- cbind(p, theta, sigma)

# B2
p <- c(0.0225, 0.255, 0.7225)
theta <- c(-4,-1, 1)
sigma <- c(0.15, 0.15, 0.15)
params <- cbind(p, theta, sigma)

# B3
p <- c(0.0225, 0.255, 0.7225)
theta <- c(-2,-0.75, 0.5)
sigma <- c(0.2, 0.2, 0.2)
params <- cbind(p, theta, sigma)

# B4
p <- c(0.0225, 0.255, 0.7225)
theta <- c(-1.5,-0.5, 0.25)
sigma <- c(0.25, 0.25, 0.25)
params <- cbind(p, theta, sigma)

# C1 - p = 0.25
p <- c(0.5625, 0.375, 0.0625)
theta <- c(-4,-0.5, 0.5)
sigma <- c(0.2, 0.1, 0.05)
params <- cbind(p, theta, sigma)

# C2
p <- c(0.5625, 0.375, 0.0625)
theta <- c(-4,-1, 1)
sigma <- c(0.3, 0.3, 0.3)
params <- cbind(p, theta, sigma)

# C3
p <- c(0.5625, 0.375, 0.0625)
theta <- c(-2,-0.75, 0.5)
sigma <- c(0.2, 0.2, 0.2)
params <- cbind(p, theta, sigma)

# C4
p <- c(0.5625, 0.375, 0.0625)
theta <- c(-1,-0.3, 0.25)
sigma <- c(0.15, 0.15, 0.15)
params <- cbind(p, theta, sigma)
