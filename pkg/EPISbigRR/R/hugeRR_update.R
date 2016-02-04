hugeRR_update <-
function(obj, update.col = c(1:length(obj$u)), Z.name, Z.index, family = gaussian(link = identity), tol.err = 1e-6, tol.conv = 1e-8, save.cache = FALSE)
{
	w <- as.numeric(obj$u^2/(1 - obj$leverage))
	w[w < tol.err] <- tol.err
	w[-update.col] = 1
	# Z.updated <- t(sqrt(w)*t(Z)) ## bug fixed 111201 -- xia
	hugeRR(y = obj$y, X = obj$X, Z.name = Z.name, Z.index = Z.index, family = family, 
		   weight = w, tol.err = tol.err, tol.conv = tol.conv, save.cache = save.cache)
}