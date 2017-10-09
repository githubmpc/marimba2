#' A function to check and adjust for label switching.
#' 
#' This function checks for label switching.
#' 
#' @param gibbs.fit Estimated Gibbs sampler output in the form of gmodel.
#' @param gparams The parameters as set in the geneticParameters function.
#' @return Estimated Gibbs sampler gmodel adjusted for label switching if required.
#' @export
#'

label.check <- function(gibbs.fit, gparams){
  
  thetas.check <- apply(gibbs.fit$chains$theta,2,mean)
  N <- gparams$N
  
  if(all(diff(thetas.check) >= 0)==T) {
      gibbs.fit <- gibbs.fit
  } else {
    index.sort <- sort.int(thetas.check, index.return=T)
    gibbs.fit$gmodel$theta <- gibbs.fit$gmodel$theta[index.sort$ix]
    gibbs.fit$gmodel$sigma <- gibbs.fit$gmodel$sigma[index.sort$ix]
    gibbs.fit$gmodel$pi <- gibbs.fit$gmodel$pi[index.sort$ix]
    gibbs.fit$gmodel$pi.child <- gibbs.fit$gmodel$pi.child[index.sort$ix]
    # need to revise following lines when make generic
    attributes(gibbs.fit$gmodel$theta)$names <- c("0","1","2")
    attributes(gibbs.fit$gmodel$sigma)$names <- c("0","1","2")
    attributes(gibbs.fit$gmodel$pi)$names <- c("0","1","2")
    attributes(gibbs.fit$gmodel$pi.child)$names <- c("0","1","2")
    
    # this bit reorders the cn
    cn.sort <- order(thetas.check) - 1
    gibbs.fit$gmodel$cn <- matrix(as.numeric(factor(as.numeric(gibbs.fit$gmodel$cn), levels=cn.sort)), nrow=N, ncol=3)
    # this line of code below will need to checked/ tested again once the model is generalised.
    gibbs.fit$gmodel$cn <-  gibbs.fit$gmodel$cn - 1
    
    y <- gibbs.fit$gmodel$y
    gmodel <- initialize_gmodel(y, gparams, 1)
    gmodel$theta <- gibbs.fit$gmodel$theta
    gmodel$sigma <- gibbs.fit$gmodel$sigma
    gmodel$pi <- gibbs.fit$gmodel$pi
    gmodel$pi.child <- gibbs.fit$gmodel$pi.child
    seed <- runif(1, min=1, max=1000000)
    set.seed(seed)
    gibbs.fit <- gibbs_genetic(gmodel, gparams)
    thetas.check.2 <- apply(gibbs.fit$chains$theta,2,mean)
    if(all(diff(thetas.check.2) >= 0)==F)  stop("this region cannot be called accurately")
  }
  return (gibbs.fit)
}
