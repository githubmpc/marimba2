#' A function to check for low ESS.
#' 
#' This function checks for low effective sample size and re-runs the gibbs sampler with a higher number of burn-ins.
#' 
#' @param gibbs.result Summarised Gibbs.sampler output obtained using gibbs.result function
#' @param gparams The parameters as set in the geneticParameters function.
#' @return Estimated Gibbs sampler gmodel adjusted for label switching if required.
#' @export
#'

ess.check <- function(gibbs.cnv){
  
  p <- gibbs.cnv$post.pi
  eff.size.pi <- as.numeric(effectiveSize(p))
  p <- gibbs.cnv$post.pi.child
  eff.size.pi.child <- as.numeric(effectiveSize(p))
  p <- gibbs.cnv$post.thetas
  eff.size.thetas <- as.numeric(effectiveSize(p))
  p <- gibbs.cnv$post.sigmas
  eff.size.sigmas <- as.numeric(effectiveSize(p))
  ess.combo <- rbind(eff.size.pi, eff.size.pi.child, eff.size.thetas, eff.size.sigmas)
  
  if(all(ess.combo >= 200)==T) {
    gibbs.cnv <- gibbs.cnv
  } else {
    gparams <- gibbs.cnv$gparams
    N <- gparams$N
    K <- gparams$K
    y <- gibbs.cnv$y
    thin <- gparams$thin
    iter <- gparams$iter
    xi <- gparams$xi
    burnin <- 1500
    model <- gparams$model
    gibbs.cnv <- gibbs.cnv.wrapper(N=N, K=K, y=y, thin=thin, iter=iter, xi=xi, burnin=burnin, model=model)
  }
  return (gibbs.cnv)
}
