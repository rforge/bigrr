bigRR.default <-
function (formula = NULL, y, X, Z, RandC = ncol(Z), data = NULL, shrink = NULL, weight = NULL,
    family = gaussian(link = identity), lambda = NULL, impute = FALSE, tol.err = 1e-6, 
	tol.conv = 1e-8, only.estimates = FALSE, GPU = FALSE, ...) 
{

    Call <- match.call()
    if (!(is.matrix(X))) 
        stop("X should be a matrix.")
    if (!(is.matrix(Z))) 
        stop("Z should be a matrix.")
    if (!(is.vector(y))) 
        stop("y should be a vector.")
	#if (GPU & !require(gputools, quietly = TRUE)) stop('Package gputools is needed for using GPU.\n')
	if (GPU) stop('GPU option is only available in R-Forge versions of "bigRR".\n')
	if (any(is.na(y))) {
		naidx <- which(is.na(y))
		y <- y[-naidx]
		X <- X[-naidx,]
		Z <- Z[-naidx,]
	}
    if (impute) {
      if (any(is.na(Z))) {
        cat('Imputing missing values...')
        nacolidx <- which(is.na(colSums(Z)))
        for (j in nacolidx) {
          naidx <- which(is.na(Z[,j]))
          Z[naidx,j] <- sample(Z[-naidx,j], length(naidx), TRUE)
        }
        cat('Done.\n')
        } else {
          cat("NOTE: no missing value exists, no need to impute.\n")
        }
    }
    N <- n <- nrow(X)
    p <- ncol(X)
    k <- ncol(Z)
    if (N != nrow(Z) | N != length(y)) 
        stop("Sizes of y, X, and Z are not all equal.")
    if (is.null(weight)) w <- rep(1, k) else w <- weight

    #G <- crossprod(sqrt(w)*t(Z)) ## bug fixed 111201 -- Xia
    wZt <- sqrt(w)*t(Z)
  ####Multiple VC code from here 131017 -- Lars
  if (length(RandC)==1) {
    if (!GPU) G <- crossprod(wZt) #else G <- gpuMatMult(t(wZt), wZt)
    ############ Bending to allow for p<n problems -- Lars (Xia added SVD)
    if (k < n) {
      eigen.values <- svd(G)$d
      min.eigen <- min(eigen.values)
      if (min.eigen < tol.err) G <- G + diag(N)*(abs(min.eigen) + tol.err) 
    }
    ############
    #invG <- solve(G)
    #L <- t(chol(G))
    svdG <- svd(G)
    L <- svdG$u %*% diag(sqrt(svdG$d))
    invG <- tcrossprod(svdG$v %*% diag(1/svdG$d), svdG$u)
    phi0 <- sa0 <- 1
    if (is.null(lambda)) {
        hm <- hglm(y = y, X = X, Z = L, family = family, conv = tol.conv) ## checked with old emme code, conv = 1e-6 removed -- Xia
    }
    else {
        start.beta = c(rep(0, p))
        start.v = c(rep(0.01, n))
        start.lambda = lambda
        start.sigma2e = 1
        cat("Only 1 iteration applied for fixed lambda")
        hm <- hglm(y = y, X = X, Z = L, family = family, startval = c(start.beta, 
                   start.v, start.lambda, start.sigma2e), maxit = 1)
    }
    phi <- as.numeric(hm$varFix)
    sa <- as.numeric(hm$varRanef)
    a <- L%*%hm$ranef
    if (!GPU) tZinvG <- crossprod(Z, invG) #else tZinvG <- gpuMatMult(t(Z), invG)
    u <- (w*tZinvG)%*%a
    qu <- GCV <- NULL
    if (!only.estimates) {
        C <- rbind(cbind(crossprod(X, X), crossprod(X, L)), cbind(crossprod(L, X), G + diag(N)*phi/sa))
        C22 <- solve(C)[(p + 1):(p + N), (p + 1):(p + N)]*phi
        if (!GPU) transf <- (w*tZinvG) %*% L #else transf <- gpuMatMult((w*tZinvG), L)
        qu <- hat.transf(C22, transf, vc = sa, w, k, N, tol.err = tol.err, GPU = GPU)
        adjsmahat <- adjbighat <- 0
        if (any(qu < tol.err) | any (qu > 1 - tol.err)) {
          adjsmahat <- sum(qu < tol.err)
          adjbighat <- sum(qu > (1 - tol.err))
        	qu[qu < tol.err] <- tol.err
        	qu[qu > (1 - tol.err)] <- 1 - tol.err
        }
        if (family$family == "gaussian") GCV <- sum(hm$resid^2)/((n - sum(hm$hv[1:n]))^2)
    }
  }
    if (!GPU & length(RandC)>1) {
      n.vc = length(RandC)
      L <- invG <- matrix(0, N, n.vc*N)
      count.L = 0
      start.col = 1
      end.col = 0
      for (n.col in RandC) {
        end.col = end.col + n.col
        count.L = count.L + 1
        G <- crossprod(wZt[start.col:end.col,])
        svdG <- svd(G)
        L[,((count.L-1)*n+1):(count.L*n)] <- svdG$u %*% diag(sqrt(svdG$d))
        invG[,((count.L-1)*n+1):(count.L*n)] <- tcrossprod(svdG$v %*% diag(1/svdG$d), svdG$u)
        start.col=end.col+1
      }
      hm <- hglm(y = y, X = X, Z = L, RandC = c(rep(N,n.vc)), family = family, conv = tol.conv, maxit=200)
      phi <- as.numeric(hm$varFix)
      sa <- as.numeric(hm$varRanef)
print(hm)
print(sa)
print(N)
      #a <- L%*%hm$ranef
      #tZinvG <- crossprod(Z, invG) 
      start.col = 1
      end.col = 0
      u <- numeric(k)
      tZinvG <- matrix(0, k, N)
      for (i.vc in 1:n.vc) {
        tmp = ((i.vc-1)*N+1):(N*i.vc)
        end.col = end.col + RandC[i.vc]
        tmp2 = start.col:end.col
        a <- L[ ,tmp] %*% hm$ranef[tmp]
        tZinvG[tmp2,] <- crossprod(Z[,tmp2], invG[,tmp]) 
        u[tmp2] <- (w[tmp2]*tZinvG[tmp2,]) %*%  a
        start.col = end.col + 1
      } 
      qu <- GCV <- NULL
      if (!only.estimates) {
        C <- rbind(cbind(crossprod(X, X), crossprod(X, L)), cbind(crossprod(L, X), crossprod(L) + diag(rep(1/sa, each = N))*phi))
        C22 <- solve(C)[(p + 1):(p + N*n.vc), (p + 1):(p + N*n.vc)]*phi
        qu = NULL
        start.col = 1
        end.col = 0
        for (i.vc in 1:n.vc) {
          tmp = ((i.vc-1)*N+1):(N*i.vc)
          end.col = end.col + RandC[i.vc]
          tmp2 = start.col:end.col
          transf <- (w[tmp2]*tZinvG[tmp2,]) %*%  L[ ,tmp]
          qu <- c(qu, hat.transf(C22[tmp,tmp], transf, vc = sa[i.vc], w = w[tmp2], k = length(tmp2), N, tol.err = tol.err, GPU = FALSE) ) 
          start.col = end.col + 1
        } 
        adjsmahat <- adjbighat <- 0
        if (any(qu < tol.err) | any (qu > 1 - tol.err)) {
          adjsmahat <- sum(qu < tol.err)
          adjbighat <- sum(qu > (1 - tol.err))
          qu[qu < tol.err] <- tol.err
          qu[qu > (1 - tol.err)] <- 1 - tol.err
        }
        #if (family$family == "gaussian") GCV <- sum(hm$resid^2)/((n - sum(hm$hv[1:n]))^2)
      }
    }
    result <- list(phi = phi, lambda = hm$varRanef, beta = hm$fixef, hglm = hm,
        u = u, leverage = qu, GCV = GCV, Call = Call, y = y, X = X, Nsmallhat = adjsmahat, Nbighat = adjbighat)
    class(result) <- "bigRR"
    return(result)
}

