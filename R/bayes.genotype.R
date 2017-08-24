#' Bayesian copy number genotyping
#' 
#' This function allows you to genotype copy number for a dataset of trios. Uses the output from the gibbsSamplerMendelian function.
#' @param model This is the output from the gibbsSamplerMendelian function.
#' @keywords marimba
#' @return The predicted copy number genotype for a dataset of trios.
#' @examples
#' classifyBayesMixture(mendelian.noerrormodel)

classifyBayesMixture <- function(model){
  cn.list <- model$classifications
  mean.0 <- Reduce("+", lapply(cn.list, function(x) x == 0)) / length(cn.list)
  mean.1 <- Reduce("+", lapply(cn.list, function(x) x == 1)) / length(cn.list)
  mean.2 <- Reduce("+", lapply(cn.list, function(x) x == 2)) / length(cn.list)
  class.bayes <- matrix(0, nrow = nrow(mean.0), ncol = ncol(mean.0))
  for( i in 1:nrow(mean.0)){
    for( j in 1:ncol(mean.0)){
      class.bayes[i,j] <- order(c(mean.0[i,j], mean.1[i,j], mean.2[i,j]), decreasing = TRUE)[1]
    }
  }
  return(class.bayes)
}