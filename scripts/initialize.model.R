#' A function to initialize the model.
#' 
#' This function allows you to initialize the simulation model with starting values for Gibbs sampling.
#' 
#' @param dat Matrix of array intensity data where columns represent mother, father and offspring.
#' @param params The parameters as set in the geneticParameters function.
#' @param comp Required default variable. Will remove in future versions of package.
#' @return Returns the set parameters
#' @export
#'

initialize_gmodel <- function(dat.y, params, comp){
  y.mf <- multi_K(dat.y, params)
  y.mf <- y.mf[[comp]]
  y.o <- dat.y[, "o"]
  y.odf <- data.frame(value = y.o)
  ## Assume that Z=0 corresponds to class 1, Z=1 corresponds to class 2
  ## Random guesses for Z
  cn.mf <- dcast(y.mf, Var1 ~ Var2, value.var="cn")[, -1]
  if(length(attributes(table(as.matrix(cn.mf))))!=3) next
  y.mf.list <- split(y.mf$value, y.mf$cn)
  ## Calculate sufficient statistics
  Ns <- sapply(y.mf.list, length)
  mus <- sapply(y.mf.list, mean, na.rm=TRUE)
  vars <- sapply(y.mf.list, var, na.rm=TRUE)
  theta <- mus
  sigma2 <- vars
  sigma <- sqrt(sigma2)
  a <- params$a
  d <- length(table(y.mf$cn))
  a <- a[1:d]
  eta1 <- params$eta1
  eta2 <- params$eta2
  pi <- rdirichlet(1, Ns+a)[1, ]
  names(pi) <- names(theta)
  tau <- rbeta(1, eta1, eta2)
  ## create initial mendelian matrix
  mendelian.probs <- gMendelian(tau.one=params$tau.one, tau.two=params$tau.two, tau.three=params$tau.three, err=params$error)
  ## Draw for offspring
  p.o <- p_offspring(cn.mf, mendelian.probs, theta)
  
  # process p.o into appropriate dimensions for < 3 components
  p.o <- dim.reduct(p.o, comp)
  
  # initialize pi.child
  pi.child <- init.pi.child(dat.y, cn.mf, params, theta)
  names(pi.child) <- names(theta)
  
  cn.prob <- cnProb(p.o, y.o, theta, sigma)
  cn.prob <- as.matrix(cn.prob)
  ##
  ## initialize copy number for offspring
  ##
  y.odf <- init.offspring.cn(cn.prob, cn.mf, comp, y.odf)
  tab <- table(y.odf$cn)
  check<-length(unique(y.mf$cn))
  if(length(tab) != check) stop("Some components unobserved")
  cn.all <- as.matrix(cbind(cn.mf, y.odf$cn))
  colnames(cn.all) <- c("m", "f", "o")
  list(y=dat.y,
       cn=cn.all,
       m.probs=mendelian.probs,
       theta=theta,
       sigma=sigma,
       pi=pi,
       pi.child=pi.child,
       tau=tau)
}
