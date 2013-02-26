bigRR <-
function(formula = NULL, y = NULL, X = NULL, Z = NULL, data = NULL, shrink = NULL, weight = NULL,
		family = gaussian(link = identity), lambda = NULL, impute = FALSE, tol.err = 1e-6, 
		tol.conv = 1e-8, only.estimates = FALSE, GPU = FALSE, ...) UseMethod("bigRR")

.onAttach <- 
		function(...)
{
	packageStartupMessage("bigRR: Fast generalized ridge regression for p >> n data")
	packageStartupMessage('Version 1.3-3 installed')
	packageStartupMessage('Authors:    Xia Shen - xia.shen@slu.se')
	packageStartupMessage('            Moudud Alam - maa@du.se')
	packageStartupMessage('            Lars Ronnegard - lrn@du.se')
	packageStartupMessage('Maintainer: Xia Shen - xia.shen@slu.se')
	packageStartupMessage('Use citation("bigRR") to know how to cite our work.')
	options(warn = -1)
}