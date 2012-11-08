core <-
function(X, y, Sigmas=NULL, ns=NULL, numdir=2, verbose=FALSE, short=TRUE,...)
{
	mf <- match.call()

	if (!is.null(Sigmas)){ p <- dim(Sigmas[[1]])[1]; h <- length(Sigmas)}

	if (is.null(Sigmas) & is.null(ns))
	{
		hlevels <- unique(y); p <- ncol(X);

		Sigmas <-list(); ns <- vector(length=length(hlevels));

		for (h in 1:length(hlevels)) 
		{
			ns[h] <- sum(y==h);
 
			tempX <- X[y==h,]; 

			Sigmas[[h]] <- t(tempX)%*%tempX/ns[h];
		}
	} 
	n <- sum(ns); ff <- ns/n; d <- numdir;

	one.core <- function(d)
	{
		if (d == 0)
		{
			Gamma.hat <- NULL;

			Sigma.hat <- matrix(0, p, p);

			for (g in 1:h) {Sigma.hat <- Sigma.hat + ff[g] * Sigmas[[g]]}
	
			term0 <- 0;

			term1 <- n/2 * log(det(Sigma.hat));

			loglik <- term0 - term1;
		}
		else if (d == p)
		{
			Gamma.hat <- diag(p);

			Sigma.hat <- matrix(0, p, p);

			for (g in 1:h){Sigma.hat <- Sigma.hat + ff[g] * Sigmas[[g]]}
	
			term0 <- 0;

			term1 <- n/2 * log(det(Sigma.hat));

			term2 <- n/2 * log(det(Sigma.hat));

			term3 <- 0;

			for (g in 1:h){term3 <- term3 + ns[g]/2 * log(det(Sigmas[[g]]))}

			loglik <- Re(term0 - term1 + term2 - term3)
		}
		else 
		{
			objfun <- function(W)
			{
				Q <- W$Qt;	d <- W$dim[1]; p <- W$dim[2];

				Sigmas <- W$Sigmas;

				n <- sum(W$ns);

				U <- matrix(Q[,1:d], ncol=d);

				V <- matrix(Q[,(d+1):p], ncol=(p-d));

				Sigma.hat <- matrix(0, p, p);

				for (g in 1:h){Sigma.hat <- Sigma.hat + ff[g] * Sigmas[[g]]}

				Ps <- projection(U, diag(p));

				# Objective function

				term0 <- 0;

				term1 <- n/2 * log(det(Sigma.hat));

				term2 <- n/2 * log(det(t(U) %*% Sigma.hat %*% U));

				term3 <- 0;

				for (g in 1:h){term3 <- term3 + ns[g]/2 * log(det(t(U) %*% Sigmas[[g]] %*% U))}

				value <- Re(term0 - term1 + term2 - term3)

				return(list(value=value))
			}
			objfun <- assign("objfun", objfun, envir=.GlobalEnv); 

			W <- list(dim=c(d, p), Sigmas=Sigmas, ns=ns);

			grassmann <- GrassmannOptim(objfun, W, verbose=verbose,...);

			Gamma.hat <- matrix(grassmann$Qt[,1:d], ncol = d);

			loglik <- tail(grassmann$fvalues, n = 1)
		}

		Sigma.hat <- matrix(0, p, p);

		for (g in 1:h){Sigma.hat <- Sigma.hat + ff[g] * Sigmas[[g]]}

		if (d != 0){Ps.hat <- projection(Gamma.hat, Sigma.hat)}

		Sigmas.hat <- list();

		for (g in 1:h)
		{
			if (d == 0){Sigmas.hat[[g]] <- Sigma.hat
			}
			else{	
				Sigmas.hat[[g]] <- Sigma.hat + t(Ps.hat) %*% ( Sigmas[[g]] - Sigma.hat) %*% Ps.hat
			}
		}
		numpar <- p*(p+1)/2 + d*(p-d) + (h-1)*d*(d+1)/2;

		aic <- -2*loglik + 2 * numpar;

		bic <- -2*loglik + log(n) * numpar;

		return(list(Gammahat=Gamma.hat, Sigmahat = Sigma.hat, Sigmashat = Sigmas.hat,
			loglik=loglik, numpar=numpar, aic=aic, bic=bic))
	}

	if (short)
	{
		fit <- one.core(d);

		ans <- list(Gammahat=fit$Gammahat, Sigmahat = fit$Sigmahat, Sigmashat = fit$Sigmashat, 
			loglik=fit$loglik, aic=fit$aic, bic=fit$bic, numpar=fit$numpar, numdir=d, 
			 model="core", call=match.call(expand.dots = TRUE), short=short)

		class(ans) <- "ldr"

		return(invisible(ans))
	}
	aic <- bic <- numpar <- loglik <- vector(length=d+1);

	Gammahat <- Sigmahat <- Sigmashat <- list();

	loglik <- numpar <- aic <- bic <- numeric(d+1);

	for (i in 0:d)
	{
		if (verbose) cat("Running CORE for d =", i, "\n")

		fit <- one.core(i)

		Gammahat[[i+1]] <-fit$Gammahat

		Sigmahat[[i+1]] <- fit$Sigmahat

		Sigmashat[[i+1]] <- fit$Sigmashat

		loglik[i+1] <- fit$loglik

		numpar[i+1] <- fit$numpar

		aic[i+1] <- fit$aic

		bic[i+1] <- fit$bic
	}
	ans <- list(Gammahat=Gammahat, Sigmahat = Sigmahat, Sigmashat = Sigmashat,
		loglik=loglik, aic=aic, bic=bic, numpar=numpar, numdir=d, model="core",
		call=match.call(expand.dots = TRUE), short=short)

	class(ans) <- "ldr"

	return(invisible(ans))
}
