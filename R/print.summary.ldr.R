print.summary.ldr <-
function(x,...)
{
  	cat("\nCall:\n"); print(x$call);

	if (x$numdir==0) return(list())

	if (identical(x$short, FALSE))
	{	
		cat("\nEstimated Basis Vectors for Central Subspace:\n");
		
		if (x$model=="pfc")
		{ 
			if ((x$structure == "aniso") | (x$structure=="unstr"))
			{
				print(round(OrthoNorm(solve(x$Deltahat[[x$numdir]])%*%x$Gammahat[[x$numdir]]), digits=4)) 

			} else {
			
				print(round(x$Gammahat[[x$numdir]], digits=4)) 
			}
		} 
		else  print(round(x$Gammahat[[x$numdir+1]], digits=4))

		cat("\nInformation Criteria:\n");

		print(x$IC);

		cat("\nLarge sample likelihood ratio test \n")

		print(x$LRT); cat("\n");

		if (x$model=="pfc")
		{
			if (!is.null(x$Rsq)) print(round(x$Rsq, digits=4)); 	cat("\n"); 
		}
	}
	else 
	{
		cat("\n\nEstimated Basis Vectors for Central Subspace:\n");

		if (x$model=="pfc")
		{ 
			if ((x$structure == "aniso") | (x$structure=="unstr"))
			{
				print(round(OrthoNorm(solve(x$Deltahat)%*%x$Gammahat), digits=4)) 
			} else{
				print(round(x$Gammahat, digits=4))
			}
		} else  print(round(x$Gammahat, digits=4))
	}
}
