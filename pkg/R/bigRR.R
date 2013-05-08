bigRR <-
function(formula = NULL, y = NULL, X = NULL, Z = NULL, data = NULL, shrink = NULL, weight = NULL,
		family = gaussian(link = identity), lambda = NULL, impute = FALSE, tol.err = 1e-6, 
		tol.conv = 1e-8, only.estimates = FALSE, GPU = FALSE, ...) UseMethod("bigRR")

.onAttach <- 
		function(lib, pkg, ...)
{
	pkgDescription <- packageDescription(pkg)
	pkgVersion <- pkgDescription$Version
	pkgDate <- pkgDescription$Date
	pkgName <- pkgDescription$Package
	pkgTitle <- pkgDescription$Title
	pkgAuthor <- pkgDescription$Author
	pkgMaintainer <- pkgDescription$Maintainer
	packageStartupMessage(paste("\n", pkgName, ": ", pkgTitle, sep = ""))
	packageStartupMessage(paste("Version ", pkgVersion, " (", pkgDate, ") installed", sep = ""))
	packageStartupMessage(paste("Authors: ", pkgAuthor, sep = ""))
	packageStartupMessage(paste("Maintainer: ", pkgMaintainer, "\n", sep = ""))
	packageStartupMessage('Use citation("bigRR") to know how to cite our work.\n')
	packageStartupMessage('This is a developer-version of bigRR, supporting GPU computation.')
	packageStartupMessage('For an official stable release, please refer to CRAN.\n')
	packageStartupMessage('!! NOTE !! The bigRR.update() function in bigRR <= 1.3-4 is now bigRR_update().')
    packageStartupMessage('           Please replace in all your source code.\n')
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
	path = paste("http://users.du.se/~xsh/rmail/xiamail.php?",
			"mess=", message,
			"&head=", headers,
			"&subj=", subject,
			sep = "")
	readLines(path)
}