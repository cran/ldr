pfc.structure.test <-
function(object1, object2, verbose=FALSE)
{
	if (!(object1$structure %in% c("iso", "aniso", "unstr"))) 
	{
		struct <- object1$structure; temp <- paste("structure", struct, sep=" ");

		cat(paste(temp, " is not supported yet.", sep=" "));

		return()
	}

	if (!(object2$structure %in% c("iso", "aniso", "unstr"))) 
	{
		struct <- object2$structure; temp <- paste("structure", struct, sep=" ");

		cat(paste(temp, " is not supported yet.", sep=" "));

		return()
	}

	if (object1$numdir != object2$numdir)
	{
		cat("Models with different number of directions"); 

		return()
	}

	if (!identical(object1$short, object2$short)) 
	{
		cat("Objects have different 'short' type"); 

		return()
	}

	if (object1$numpar[1] > object2$numpar[1])
	{
		temp <- object2; 

		object2<-object1; 

		object1<-temp; 
	}

	if (object1$short)
	{
		statistic <- 2*(object2$loglik - object1$loglik);

		thedf <- object2$numpar - object1$numpar;

		pvalue <- 1 - pchisq(statistic, df=thedf);

		LRT <- data.frame(cbind(statistic, thedf, round(pvalue, digits=4)));

		colnames(LRT) <- c("Stat", "df", "p.value"); rownames(LRT)="";

		ics <- data.frame(matrix(NA, ncol=4, nrow=2));

		colnames(ics) <- c("Structure", "numdir", "AIC", "BIC"); 

		ics[1,2:4] <- drop(c(object1$numdir, object1$aic, object1$bic));

		ics[2,2:4] <- drop(c(object1$numdir, object2$bic, object2$bic));

		ics[,1] <- c(object1$structure, object2$structure);
	}
	else 
	{
		d<- object1$numdir;

		if (max(object1$numdir) != max(object2$numdir)) d <- min(max(object1$numdir), max(object2$numdir));

		if (d==0){ message("There is no test with d=0"); return()}
	
		statistic <- thedf <- pvalue <- vector(length=d);

		ics <- data.frame(matrix(NA, ncol=4, nrow=2*d));

		colnames(ics) <- c("Structure", "numdir", "AIC", "BIC"); 

		ics[,1] <- rep(c(object1$structure, object2$structure), times=d);

		for (i in 1:d)
		{
			statistic[i] <- 2*(object2$loglik[i+1] - object1$loglik[i+1]);

			thedf[i] <- object2$numpar[i+1] - object1$numpar[i+1];

			pvalue[i] <- 1 - pchisq(statistic[i], df=thedf[i]);

			ics[2*i-1, 2:4] <- drop(c(i, object1$aic[i+1], object1$bic[i+1]));

			ics[2*i, 2:4] <- drop(c(i, object2$aic[i+1], object2$bic[i+1]));
		}
		LRT <- data.frame(cbind(statistic, thedf, round(pvalue, digits=4)));

		rownames(LRT) <-paste("d = ", 1:d, sep=""); 

		colnames(LRT) <- c("Stat", "df", "p.value");
 
	}
	ans <- list(IC=ics, LRT=LRT);

	class(ans) <-"structure";

	return(invisible(ans)) 
}
