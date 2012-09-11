OrthoNorm <-
function(X)
{
	X <- as.matrix(X); n<-nrow(X); p<-ncol(X); M <- NULL;

	M <- cbind(M, X[,1]);

	if(p > 1)
	{
     		for(k in 2:p) 
		{
			one <- rep(0, n);

			for(i in 1:(k-1)) 
			{
				oneki <- as.vector((t(M[,i]) %*% X[,k])/(t(M[,i]) %*% M[,i]));

				one <- one + oneki * M[,i];
			}
			M <- cbind(M, X[,k] - one);
		}
		M <- round(M, digits=4);
	} 
	return(apply(M, 2, function(x) x/sqrt(t(x)%*%x)))
}
