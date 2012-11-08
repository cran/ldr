print.ldr <-
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
	}
	else 
	{
		cat("\n\nEstimated Basis Vectors for Central Subspace:\n");

        	if (x$model=="pfc")
		{ 
			if ((x$structure == "aniso") | (x$structure=="unstr"))
			{
				print(round(OrthoNorm(solve(x$Deltahat)%*%x$Gammahat), digits=4)) 
			} else {
				print(round(x$Gammahat, digits=4)); cat("\n"); 
			}
		} else print(round(x$Gammahat, digits=4)); cat("\n"); 
	}
}
