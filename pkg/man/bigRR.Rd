\name{bigRR}
\alias{bigRR}
\alias{bigRR.default}
\alias{bigRR.formula}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Fitting big ridge regression}
\description{
  Function fits big ridge regression with special computational advantage for the cases when number of shrinkage parameter exceeds number of observations. The shrinkage parameter, lambda, can be pre-specified or estimated along with the model. Any subset of model parameter can be shrunked.   

}
\usage{
bigRR(formula = NULL, y = NULL, X = NULL, Z = NULL, data = NULL, shrink = NULL, lambda = NULL, tol.err = 1e-6, only.estimates = FALSE, ...)
\method{bigRR}{default}(formula = NULL, y , X , Z , data = NULL, shrink = NULL, lambda = NULL, tol.err = 1e-6, only.estimates = FALSE, ...)
\method{bigRR}{formula}(formula, y = NULL, X = NULL, Z = NULL, data = NULL, shrink = NULL, lambda = NULL, tol.err = 1e-6, only.estimates = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ a two-sided model \code{formula}}
  \item{y}{ response variable; either y or \code{formula} is to be used.}
  \item{X}{ design matrix related to the parameters not to be shrunk; not redundant if formula is already used }
  \item{Z}{ design matrix associated with shrinkage parameters to be shrunk; not required if model formula is used}
  \item{data}{ a data frame}
  \item{shrink}{ either a numeric or a character vector specifying the position or the names if the variables whose coefficients are to be shrunken.}
  \item{lambda}{ the shrinkage parameter determines the amount of shrinkage. Default is NULL meaning that it is to be estimated along with other model parameters. }
  \item{tol.err}{ tolerance level; default value is 1e-8 }
  \item{only.estimates}{logical; TRUE if hat values are to be returned. Default is FALSE and the hat values are not computed.}
  \item{\dots}{ unused arguments }
}
\details{
  It fits ridge regression using the algorithm presented in Ronnegard et al. (2010). The computationla intensity of the estimation depends on the number of parameters but not the number of observatons. 
The model can be specified by using \code{formula} or by design matrices. If both are used then the formula interface will be used. The shrink argument specifies the subset of parameters  to be estimated. If the model is specified by using formula and a dot (.) is used in the right hand side of the formula then \code{shrink} indicates the variables in the data frame. Otherwise it represents the respective variable in the model formula. It is ignored it model is specified by using design matrices.  
}
\value{
  Returns a lis of object class \code{bigRR} containing the following values 
 
  \item{phi }{error variance}
  \item{lambda }{usual \code{labda} parameter of the gidre regression representing the amount of shrinkage}
  \item{beta}{subset of model parameters which is/are not shrunken, i.e. those associated with X matrix}
  \item{u}{subset of model parameters which are shrunken, i.e. those associated with Z matrix}
  \item{leverage}{hat values}
  \item{Call}{how the bigRR was called.}
}
\references{
Shen, X., Alam, M., Fikse, F. and Ronnegard, L. (2012). \emph{Fast generalized ridge regression for models including heteroscedastic effects in quantitative genetics}. Manuscript.
}
\author{Xia Shen, Moudud Alam, Lars Ronnegard}

\seealso{ \code{lm.ridge} in MASS library.}

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

