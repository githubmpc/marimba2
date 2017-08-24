#' Simulating copy number for a collection of trios.
#' 
#' This function allows you to simulate copy number for a dataset of trios. This can be used as input for the gibbsSamplerMendelian function.
#' @param n The number of trios.
#' @param p A vector representing the probabilities of each copy number state.
#' @param theta A vector representing the mean LRR of each copy number state.
#' @param sigma A vector representing the average LRR variance for each copy number state.
#' @keywords marimba
#' @return A collection of trios with their copy number genotypes. 
#' @examples
#' data.1 <- simulate_data(66, p = p.1000gp, theta = theta.1000gp, sigma = sigma.1000gp)
#' saveRDS(data.1, "simulatedData.rds")

cn.simulate.data <- function(n, p, theta, sigma){
  # set seed so we get the same results for given parameter values
  set.seed(668899)
  c_m <- sample(1:3, size = n, replace = TRUE, prob = p)
  c_f <- sample(1:3, size = n, replace = TRUE, prob = p)
  c_o <- rep(0, length = n)
  for(i in 1:n){
    cn_m <- c_m[i]
    cn_f <- c_f[i] 
    p.offspring <- mendelian.probs[,cn_m, cn_f]
    c_o[i] <- sample(1:3, size = 1, prob = p.offspring)
  }
  y_m <- rnorm(n, mean = theta[c_m], sd = sigma[c_m])
  y_f <- rnorm(n, mean = theta[c_f], sd = sigma[c_f])
  y_o <- rnorm(n, mean = theta[c_o], sd = sigma[c_o])
  
  y.mat <- cbind(cbind(y_m, y_f), y_o)
  cn.mat <- cbind(cbind(c_m, c_f), c_o)
  return(list(response = y.mat, cn = cn.mat))
}