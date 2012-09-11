summary.ldr <-
function(object, y=NULL, ...)
{
	mf <- match.call();

	lrtest <- function(object)
	{
		Tests <- Dfs <- Pvals <-cnames <- vector(length=object$numdir);

		for (w in 1:object$numdir)
		{
			Tests[w] <- 2*(object$loglik[object$numdir+1]-object$loglik[w]);

			Dfs[w] <- object$numpar[object$numdir+1]-object$numpar[w];

			Pvals[w] <- 1-pchisq(Tests[w], df=Dfs[w]);

			cnames[w] <- paste(paste(w-1, "D vs >= ", sep=""), paste(w, "D", sep=""), sep="");
		}
		ans <- data.frame(cbind(Tests, Dfs, Pvals)); rownames(ans) <- cnames;

		colnames(ans) <- c("Stat", "df", "p.value");

		return(ans)
	}

	ics <- function(object)
	{
		mat <- data.frame(rbind(object$aic, object$bic));

		v1 <- c("aic", "bic"); 

		v2 <-paste("d=", 0:(object$numdir), sep="");

		dimnames(mat) = list(v1, v2);

		return(mat)
	}

	r2 <- function(object, y)
	{
		R2 <- cnames <- vector(length=object$numdir); 

		for (w in 1:object$numdir)
		{
			if (identical(object$short, TRUE))
			{
				tempd <-data.frame(cbind(y-mean(y), object$Xc%*%solve(object$Deltahat)%*%matrix(object$Gammahat[, 1:w], ncol=w))); 

				colnames(tempd)[1]<-"shorty";
			}
			if (identical(object$short, FALSE))
			{
				tempd <-data.frame(cbind(y-mean(y), object$Xc%*%solve(object$Deltahat[[w]])%*%object$Gammahat[[w]])); 

				colnames(tempd)[1]<-"shorty";
			}
			R2[w] <- summary(lm(shorty~.-1, data=tempd))$r.squared;
		}
		Rsq <- data.frame(rbind(object$evalues, R2));

		v1 <- c("Eigenvalues", "R^2(OLS|pfc)"); 

		v2 <-paste("Dir", 1:object$numdir, sep="");

		dimnames(Rsq) = list(v1, v2);

		return(Rsq)
	}

	ans <- object; class(ans) <-"summary.ldr";

	if (identical(object$short, FALSE))
	{	
		ans$LRT <- lrtest(object);

		ans$IC <- ics(object);

		if (object$model=="pfc") 
		{
			if (is.null(y)) ans$Rsq <- NULL else ans$Rsq <- r2(object, y=as.numeric(y));
		}

		return(ans);
	}
	return(ans)
}
