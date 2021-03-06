% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{initialize_chains}
\alias{initialize_chains}
\title{A Gibbs Sampler for genotyping copy number in trios}
\usage{
initialize_chains(states, N, K, S)
}
\arguments{
\item{S}{The number of iterations for the Gibbs Sampler (suggest default 10000).}

\item{y}{Data$response where data is as generated by the simulate_data function.}

\item{mu}{Defaults to -.02.}

\item{xi}{Defaults to sqrt(200).}

\item{nu}{Defaults to 3.}

\item{lambda}{Defaults to 0.05.}

\item{alpha}{Alpha1 and alpha2 and psi refer to HWE parameters as specified in Cardin et al (alpha=w, psi=lambda).}
}
\description{
This function allows you to obtain a posterior probability distribution for copy number genotyping on Illumina arrays for an fixed genetic model where tau1=1, tau2=0.5 and tau3=0
}
\keyword{marimba}
