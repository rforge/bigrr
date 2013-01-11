bigRR <-
function(formula = NULL, y = NULL, X = NULL, Z = NULL, data = NULL, shrink = NULL,
lambda = NULL, tol.err = 1e-6, only.estimates = FALSE, ...) UseMethod("bigRR")

.onAttach <- 
		function(...)
{
	cat("bigRR: Fast generalized ridge regression for p >> n data\n")
	cat('Version 1.3-3 installed\n')
	cat('\n')
	cat('Authors:    Xia Shen - xia.shen@slu.se\n')
	cat('            Moudud Alam - maa@du.se\n')
	cat('            Lars Ronnegard - lrn@du.se\n')
	cat('\n')
	cat('Maintainer: Xia Shen - xia.shen@slu.se\n')
	options(warn = -1)
}