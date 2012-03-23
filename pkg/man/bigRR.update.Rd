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
\examples{
# --------------------------------------------- #  
# Arabidopsis example #
# --------------------------------------------- #  
\dontrun{
require(bigRR)
data(Arabidopsis)
X <- matrix(1, length(y), 1)
GBLUP.result <- bigRR(y = y, X = X, Z = scale(Z), family = binomial(link = 'logit'))
HEM.result <- bigRR.update(GBLUP.result, scale(Z))
split.screen(c(1, 2))
split.screen(c(2, 1), screen = 1)
screen(3); plot(abs(GBLUP.result$u), cex = .6, col = 'slateblue')
screen(4); plot(abs(HEM.result$u), cex = .6, col = 'olivedrab')
screen(2); plot(abs(GBLUP.result$u), abs(HEM.result$u), cex = .6, pch = 19, col = 'darkmagenta')
}
}
