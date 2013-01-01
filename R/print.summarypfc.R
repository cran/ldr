print.summarypfc <-
function(x,...)
{
	orthonorm <- function(u) 
	{
		if (is.null(u)) return(NULL)
		if (!(is.matrix(u))) u <- as.matrix(u);
		dd <- dim(u); n <- dd[1]; p <-dd[2];

		if (prod(abs(La.svd(u)$d) > 1e-08) == 0) stop("collinears vectors in orthonorm")
		if (n < p)
		{
			warning("There are too much vectors to orthonormalize in orthonorm.")
			u <- as.matrix(u[, 1:p])
			n <- p
    		}
    		v <- u;
    		if (p > 1)
		{
			for (i in 2:p)
			{
				coef.proj <- c(crossprod(u[, i], v[, 1:(i - 1)]))/diag(crossprod(v[, 1:(i - 1)]));
				v[, i] <- u[, i] - matrix(v[, 1:(i - 1)], nrow = n) %*% matrix(coef.proj, nrow = i - 1)
			}
    		}
		coef.proj <- 1/sqrt(diag(crossprod(v)))
		return(t(t(v) * coef.proj))
	}
  	cat("\nCall:\n"); print(x$call);

	ans <- summary.pfc(x);

	if (identical(x$numdir.test, TRUE))
	{	
		cat("\nEstimated Basis Vectors for Central Subspace:\n");

		if ((x$structure == "aniso") | (x$structure=="unstr"))
		{
			print(round(orthonorm(solve(x$Deltahat[[x$numdir]])%*%x$Gammahat[[x$numdir]]), digits=4)) 
		} else {
			
			print(round(x$Gammahat[[x$numdir]], digits=4)) 
		}
		
		cat("\nInformation Criterion:\n");
		print(ans$IC);

		cat("\nLarge sample likelihood ratio test \n")
		print(ans$LRT); cat("\n");

		if (is.numeric(x$y)) print(round(ans$Rsq, digits=4)); 	cat("\n");  	
	}
	else 
	{
		cat("\n\nEstimated Basis Vectors for Central Subspace:\n");

		if ((x$structure == "aniso") | (x$structure=="unstr"))
		{
			print(round(orthonorm(solve(x$Deltahat)%*%x$Gammahat), digits=4)) 
		} else {
			
			print(round(x$Gammahat, digits=4)); cat("\n"); 
		}
	}
}
