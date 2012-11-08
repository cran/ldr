ldr <-
function(X, y=NULL, fy=NULL, Sigmas=NULL, ns=NULL, ycat=TRUE, numdir=2, model=c("core", "lad", "pfc"), structure="iso", verbose=FALSE, short=TRUE,...)
{
	if (model=="pfc")
	{	
		if (is.null(fy)){stop("fy is not provided"); return()}
 
 		return(invisible(pfc(X=X, fy=fy, numdir=numdir, structure=structure, verbose=verbose, short=short, ...)))
	}

	if (model=="lad") 
	{
		if (is.null(y)){stop("The response is needed"); return()}
	
		return(invisible(lad(X=X, y=y, ycat=ycat, numdir=numdir, verbose=verbose, short=short,...)))
	}

	if (model=="core") 
	{
		if (!is.null(Sigmas) & !is.null(ns)) return(invisible(core(Sigmas=Sigmas, ns=ns, numdir=numdir, verbose=verbose, short=short,...)))
        
        	return(invisible(core(X=X, y=y, numdir=numdir, verbose=verbose, short=short,...)))
	}
}
