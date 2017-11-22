#' Simulating copy number for a collection of trios.
#' @param N number of trios to simulate
#' @param error probability of non-Mendelian event in offspring
#' @param params A dataframe of the parameters of pi, thetas and sigmas.
#' @keywords marimba2
#' @return A collection of trios with their copy number genotypes. 
#' @export

.simulate_data_multi <- function(params, n, mendelian.probs, GP){
  gp <- GP
  states <- gp$states + 1
  K <- nrow(params)
  p <- params$p
  theta <- params$theta
  sigma <- params$sigma
  c_m <- sample(1:K, size = n, replace = TRUE, prob = p)
  c_f <- sample(1:K, size = n, replace = TRUE, prob = p)
  c_o <- rep(NA, length = n)
  M <- cn_adjust2(gp)
  c_m <- c_m + M
  c_f <- c_f + M
  for(i in 1:n){
    cn_m <- c_m[i] + 1
    cn_f <- c_f[i] + 1
    p.offspring <- mendelian.probs[, cn_m, cn_f]
    p.offspring <- p.offspring[states]
    c_o[i] <- sample(1:K, size = 1, prob = p.offspring)
  }
  c_o <- c_o + M
  id.index <- formatC(seq_len(n), flag="0", width=3)
  logr.tbl <- tibble(m=rnorm(n, mean = theta[c_m], sd = sigma[c_m]),
                     f=rnorm(n, mean = theta[c_f], sd = sigma[c_f]),
                     o=rnorm(n, mean = theta[c_o], sd = sigma[c_o]),
                     id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="log_ratio", -id) 
  cn.mat <- cbind(c_m, c_f, c_o)
  colnames(cn.mat) <- c("m", "f", "o")
  cn.tbl <- as.tibble(cn.mat) %>%
    mutate(id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="copy_number", -id)
  tbl <- left_join(logr.tbl, cn.tbl, by=c("id", "family_member")) %>%
    mutate(family_member=factor(family_member, levels=c("f", "m", "o"))) %>%
    arrange(id, family_member)
  tbl
}

simulate_data_multi <- function(params, N, error=0, GP){
  gp <- GP
  mendelian.probs <- gMendelian.multi()
  tbl <- .simulate_data_multi(params, N, mendelian.probs, GP)
  ## standardize
  tbl <- tbl %>%
    mutate(log_ratio=(log_ratio-median(log_ratio))/sd(log_ratio))
  M <- cn_adjust2(gp)
  z <- tbl$copy_number - M
  ##
  ## update parameters to be same as empirical values
  ##
  stats <- component_stats(tbl)
  params$p <- stats$n/sum(stats$n)
  params$sigma <- stats$sd
  params$theta <- stats$mean
  loglik <- sum(dnorm(tbl$log_ratio, params$theta[z],
                      params$sigma[z], log=TRUE))
  truth <- list(data=tbl, params=params, loglik=loglik)
  truth
}

.simulate_data <- function(params, n, mendelian.probs){
  p <- params$p
  theta <- params$theta
  sigma <- params$sigma
  c_m <- sample(1:3, size = n, replace = TRUE, prob = p)
  c_f <- sample(1:3, size = n, replace = TRUE, prob = p)
  c_o <- rep(NA, length = n)
  for(i in 1:n){
    cn_m <- c_m[i]
    cn_f <- c_f[i]
    p.offspring <- mendelian.probs[, cn_m, cn_f]
    c_o[i] <- sample(1:3, size = 1, prob = p.offspring)
  }
  id.index <- formatC(seq_len(n), flag="0", width=3)
  logr.tbl <- tibble(m=rnorm(n, mean = theta[c_m], sd = sigma[c_m]),
                     f=rnorm(n, mean = theta[c_f], sd = sigma[c_f]),
                     o=rnorm(n, mean = theta[c_o], sd = sigma[c_o]),
                     id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="log_ratio", -id) 
  cn.mat <- cbind(c_m, c_f, c_o)
  colnames(cn.mat) <- c("m", "f", "o")
  cn.mat <- cn.mat - 1
  cn.tbl <- as.tibble(cn.mat) %>%
    mutate(id=factor(paste0("trio_", id.index))) %>%
    gather(key="family_member", value="copy_number", -id)
  tbl <- left_join(logr.tbl, cn.tbl, by=c("id", "family_member")) %>%
    mutate(family_member=factor(family_member, levels=c("f", "m", "o"))) %>%
    arrange(id, family_member)
  tbl
}

simulate_data <- function(params, N, error=0){
  mendelian.probs <- gMendelian()
  tbl <- .simulate_data(params, N, mendelian.probs)
  ## standardize
  tbl <- tbl %>%
    mutate(log_ratio=(log_ratio-median(log_ratio))/sd(log_ratio))
  z <- tbl$copy_number + 1
  ##
  ## update parameters to be same as empirical values
  ##
  stats <- component_stats(tbl)
  params$p <- stats$n/sum(stats$n)
  params$sigma <- stats$sd
  params$theta <- stats$mean
  loglik <- sum(dnorm(tbl$log_ratio, params$theta[z],
                      params$sigma[z], log=TRUE))
  truth <- list(data=tbl, params=params, loglik=loglik)
  truth
}