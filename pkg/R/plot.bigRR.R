plot.bigRR <-
function(obj, alpha = .36, ...) {
	par(mfrow = c(3, 1))
	plot(obj$u, main = 'Shrinkage Estimates', ylab = 'Effect', xlab = 'Index', 
		 col = rgb(1, 0, 0, alpha), type = 'h')
	plot(abs(obj$u), main = paste('Intra-class Correlation:', round(obj$lambda/(obj$lambda + obj$phi), digits = 6)), 
	     ylab = '|Effect|', xlab = 'Index', col = rgb(1, 0, 1, alpha), type = 'h')
 	plot(obj$leverage, main = 'Leverages', ylab = 'Size', xlab = 'Index', col = rgb(0, 0, 1, alpha), pch = 19, cex = .36)
}

