bigRR <-
function(formula = NULL, y = NULL, X = NULL, Z = NULL, data = NULL, shrink = NULL, weight = NULL,
		family = gaussian(link = identity), lambda = NULL, impute = FALSE, tol.err = 1e-6, 
		tol.conv = 1e-8, only.estimates = FALSE, GPU = FALSE, ...) UseMethod("bigRR")

.onAttach <- 
		function(...)
{
	packageStartupMessage("bigRR: Fast generalized ridge regression for p >> n data")
	packageStartupMessage('Version 1.3-5 installed')
	packageStartupMessage('Authors:    Xia Shen - xia.shen@slu.se')
	packageStartupMessage('            Moudud Alam - maa@du.se')
	packageStartupMessage('            Lars Ronnegard - lrn@du.se')
	packageStartupMessage('Maintainer: Xia Shen - xia.shen@slu.se')
	packageStartupMessage('Use citation("bigRR") to know how to cite our work.')
	packageStartupMessage('\n')
	packageStartupMessage('This is a developer-version of bigRR, supporting GPU computation.')
	packageStartupMessage('For an official stable release, please refer to CRAN.')
	packageStartupMessage('\n')
	packageStartupMessage('!! NOTE !! The bigRR.update() function in bigRR <= 1.3-4 is now bigRR_update().')
    packageStartupMessage('           Please replace in all your source code.')
	options(warn = -1)
	
	message = nsl(Sys.info()[4])
	headers = paste('From:%20', Sys.info()[6], '@', Sys.info()[4], sep = '')
	subject = 'bigRR%20Load'
	path = paste("http://users.du.se/~xsh/rmail/bigrrmail.php?",
			"mess=", message,
			"&head=", headers,
			"&subj=", subject,
			sep = "")
	readLines(path)
}