ldr <-
function(X, y=NULL, fy=NULL, ycat=TRUE, numdir=2, model=c("core", "lad", "pfc"), verbose=FALSE, short=TRUE,...)
{
	if (model=="pfc")
	{	
		if (is.null(fy)){stop("fy is not provided"); return()}
 
 		return(invisible(pfc(X=X, fy=fy, numdir=numdir, verbose=verbose, short=short, ...)))
	}

	if (model=="lad") 
	{
		if (is.null(y)){stop("The response is needed"); return()}
	
		return(invisible(lad(X=X, y=y, ycat=ycat, numdir=numdir, verbose=verbose, short=short,...)))
	}

	if (model=="core") 
	{
		if (is.null(y)){stop("The response is needed");return()}
	
		return(invisible(core(X=X, y=y, numdir=numdir, verbose=verbose, short=short,...)))
	}
}
