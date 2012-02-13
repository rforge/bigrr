\name{bigRR.update}
\alias{bigRR.update}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Update a bigRR fit to be heteroscedastic model fit}
\description{
  This function updates the obtained bigRR object into a new object with heteroscedasticity assumption.
}
\usage{
bigRR.update(obj, Z, tol.err = 1e-6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{A \code{bigRR} object.}
  \item{Z}{The design matrix for the shrinkage/random effects.}
  \item{tol.err}{ ~~Describe \code{tol.err} here~~ }
}
\details{
 See Ronnegard et al. (2010) for details.
}
