print.bigRR <-
function(x, ...) {
	Obj <- x
	opt <- options()
	dgt <- opt$digits
	cat("Call: \n")
	print(Obj$Call)
	cat("\n")
	cat("Non-shrinkage parameter estimates:\n")
	print(Obj$beta)
	cat("\n")
	cat("Shrinkage parameter estimates:\n")
	print(Obj$u)
	print(paste("lambda =", round(Obj$lambda, dgt)))
}

