screen.pfc <-
function(X, fy)
{
	unipfc = function(x)
	{
		temp.dat<-data.frame(cbind(x, fy)); 

		xnam<-paste("xx", 1:r, sep="");

		names(temp.dat)<-c("yy", xnam);

		fm.lm<- as.formula( paste("yy ~ ", paste(xnam, collapse= "+")));

		fstat<-summary( lm(fm.lm, data=temp.dat))$fstatistic;

		return( as.vector( c(fstat[1], 1-pf(fstat[1], df1=fstat[2], df2=fstat[3]))))
	}
	nobs <- dim(X)[1]; 

	npred <- dim(X)[2]; 

	r <- dim(fy)[2]; 

	if (is.null(colnames(X))) colnames(X) <- paste(X, 1:npred, sep="")

	cnames <- c("F", "Pvalue", "Index");

	OneXs<-data.frame(matrix(data=NA, ncol=3, nrow=npred)); 

	colnames(OneXs) <- cnames;

	rownames(OneXs) <- colnames(X);

	for(i in 1:npred)
	{	
		tempX<-array(X[,i], c(nobs, 1)); 

		OneXs[i,]<-c(round(unipfc(x=tempX), digits=5), i);
	}
	o <- order(as.numeric(OneXs$Pvalue), decreasing=FALSE)

	return(OneXs[o,]);
}
