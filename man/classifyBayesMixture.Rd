% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bayes.genotype.R
\name{classifyBayesMixture}
\alias{classifyBayesMixture}
\title{Bayesian copy number genotyping}
\usage{
classifyBayesMixture(model)
}
\arguments{
\item{model}{This is the output from the gibbsSamplerMendelian function.}
}
\value{
The predicted copy number genotype for a dataset of trios.
}
\description{
This function allows you to genotype copy number for a dataset of trios. Uses the output from the gibbsSamplerMendelian function.
}
\examples{
classifyBayesMixture(mendelian.noerrormodel)
}
\keyword{marimba}
