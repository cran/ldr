lad <-
function(X, y, ycat=TRUE, numdir=2, nslices=NULL, verbose=FALSE, short=TRUE,...)
{
	#Initialization
	mf <- match.call();

	op <- dim(X); p <- op[2]; n <- op[1];
	
	if (p > n) stop("This method requires n larger than p");

	Xc <- scale(X, TRUE , FALSE); 

	Sigmatilde <- t(Xc)%*%Xc/n;

	d <- numdir;

	if (ycat)
	{	groups <- unique(sort(y)); 

		nslices <- length(groups); 

		N_y <-vector(length=nslices);
 
		Deltatilde_y = vector(length=nslices, "list");

		for (i in 1:nslices)
		{	
			N_y[i] <- length(y[y==groups[i]]);

			X_y <- X[which(y == groups[i]),]; 

			X_yc <- scale(X_y, TRUE, FALSE) ; 

			Deltatilde_y[[i]] <- t(X_yc)%*%(X_yc)/N_y[i];
		}		
	} 
	if (!ycat) 
	{	
		ysort <- ldr.slices(y, nslices);

		N_y <-vector(length=nslices); 

		Deltatilde_y <- vector(length=nslices, "list");

		for (i in 1:nslices)
		{	
			N_y[i] <- sum(ysort$slice.indicator==i);

			X_y <- X[which(ysort$slice.indicator == i),]; 

			X_yc <- scale(X_y, TRUE, FALSE); 

			Deltatilde_y[[i]] <- t(X_yc)%*%(X_yc)/N_y[i];
		}
	}
	Sigmas <- list(Sigmatilde=Sigmatilde, Deltatilde_y=Deltatilde_y, N_y=N_y);

  	Deltatilde <- 0;

	for (i in 1:nslices) Deltatilde <- Deltatilde + N_y[i] * Deltatilde_y[[i]];

	onelad <-function(d)
	{
		if (d==0)
		{
			Gammahat <- NULL;

			Deltahat <- Deltatilde;

			terme0 <- -(n*p)*(1+log(2*pi))/2;

			terme1 <- -n*log( det(Sigmatilde) )/2; 

			loglik = terme0 + terme1;

			Deltahat_y <- Deltatilde_y;

			numpar <- p + p*(p+1)/2;

			AIC <- -2*loglik + 2 * numpar;

			BIC <- -2*loglik + log(n) * numpar;	
		}
		else if (d==p)
		{
			Gammahat=diag(p); 

			Deltahat=Deltatilde;

			terme0 <- -(n*p)*(1+log(2*pi))/2;

			terme3<-0;

			for (i in 1:length(N_y))
			{	
				terme3<- terme3 - N_y[i] * log( det(Deltatilde_y[[i]]) )/2;
			} 
			loglik = terme0 + terme3;

			Deltahat_y <- Deltatilde_y;

			numpar <- p+(nslices-1)*d + p*(p+1)/2 + (nslices-1)*d*(d+1)/2; 

			AIC <- -2*loglik + 2 * numpar;

			BIC <- -2*loglik + log(n) * numpar;
		}
		else 
		{
			objfun <- function(W)
			{	
				Q <- W$Qt; d <- W$dim[1]; p <- W$dim[2]; 

				Sigmas <- W$Sigmas; 

				n <- sum(Sigmas$N_y); 
 
				U <- matrix(Q[,1:d], ncol=d);

				V <- matrix(Q[,(d+1):p], ncol=(p-d));

				# Objective function

				terme0 <- -(n*p)*(1+log(2*pi))/2;

				terme1 <- -n*log(det(Sigmas$Sigmatilde))/2; 

				term <- det(t(U)%*%Sigmas$Sigmatilde%*%U);

				if (is.null(term)|is.na(term)|(term <= 0)) return(NA) else terme2 <- n*log(term)/2 ;

				terme3<-0;

				for (i in 1:length(Sigmas$N_y))
				{
					term <- det(t(U)%*%Sigmas$Deltatilde_y[[i]]%*%U); 

					if (is.null(term)|is.na(term)|(term <=0)) {break; return(NA)};

					terme3 <- terme3 - Sigmas$N_y[i] * log(term)/2;
				} 
				value <- terme0 + terme1 + terme2 + terme3;
		
				#Gradient 

				terme0 <- solve(t(U)%*%Sigmas$Sigmatilde%*%U);

				terme1 <- sum(Sigmas$N_y) * t(V)%*%Sigmas$Sigmatilde%*%U%*%terme0;

				terme2 <- 0;

				for (i in 1:length(Sigmas$N_y))
				{
					temp0 <- solve(t(U)%*%Sigmas$Deltatilde_y[[i]]%*%U);

					temp <-Sigmas$N_y[i]*t(V)%*%Sigmas$Deltatilde_y[[i]]%*%U%*%temp0;

					terme2 <- terme2 - temp;
				}
				gradient <- t(terme1 + terme2);

				return(list(value=value, gradient=gradient))
			}
			objfun <- assign("objfun", objfun, envir=.GlobalEnv); 

			InitialMatrix <- function(X, y, ycat, nslices)
			{
				op <- dim(X); n <- op[1]; p <- op[2];

				if (ycat) 
				{
					yunique <- sort(unique(y)); 

					ygroups <- y; 

					nslices=length(unique(y)); 
				}
				else
				{
					Slicing <- ldr.slices(y, nslices=nslices);

					yunique <- sort(unique(Slicing$slice.indicator));

					ygroups <- Slicing$slice.indicator;
				}
				SigmaX <- var(X);

				X0 <- as.matrix(X - kronecker(cbind(rep(1,n)),rbind(apply(X,2,mean))));

				svdS<-svd(SigmaX); 

				InvSqrtX  <- svdS$u %*% diag(c(svdS$d)^(-0.5)) %*% t(svdS$v);

				Z <- X0%*%InvSqrtX;

				pv <- vector(length=nslices);

				for (i in 1:nslices) pv[i] <- sum(ygroups==yunique[i])/length(y);

				cez <- matrix(0,p,nslices);

				cezz <- array(0,dim=c(p,p,nslices));

				for (i in 1:nslices)
				{	
					for (j in 1:n)
					{	
						cez[,i] <- cez[,i] + Z[j,]*(ygroups[j]==i);

        					cezz[,,i] <- cezz[,,i] + Z[j,]%*%t(Z[j,])*(ygroups[j]==i);
    					}
    					cez[,i] <- cez[,i]/n/pv[i];

    					cezz[,,i] <- cezz[,,i]/n/pv[i];
				}
				mat1 <- mat2 <- mat4 <- mat5 <- matrix(0,p,p); mat3 <- 0;  

				for (i in 1:nslices)
				{	
					mat1 <- mat1 + cezz[,,i]%*%cezz[,,i]*pv[i];

					mat2 <- mat2 + cez[,i]%*%t(cez[,i])*pv[i];

					mat3 <- mat3 + sum(cez[,i]^2)*pv[i];

					mat4 <- mat4 + cezz[,,i]%*%cez[,i]%*%t(cez[,i])*pv[i];

					mat5 <- mat5 + sum(cez[,i]^2)*cez[,i]%*%t(cez[,i])*pv[i];
				}
				# Here are the directions
				# Using DR

				dr <- 2*mat1 + 2*mat2%*%mat2 + 2*mat3*mat2 - 2*diag(p);

				vector.dr <- eigen(dr)$vectors;

				dr.direction <- InvSqrtX %*%vector.dr;

				Wdr <- OrthoNorm(dr.direction);

				# using SIR
				sir <- mat2;  vector.sir<-eigen(sir)$vectors;

				sir.direction<-InvSqrtX %*%vector.sir;

				Wsir <- OrthoNorm(sir.direction);

				# Using SAVE
				save <- mat1 - diag(p) - mat4 - t(mat4) + 2*sir + mat5;

				vector.save <- eigen(save)$vectors;

				save.direction <- InvSqrtX %*%vector.save;

				Wsave <- OrthoNorm(save.direction);

				return(list(Wsir=Wsir, Wsave=Wsave, Wdr=Wdr))
			}

			if (is.null(mf$sim_anneal))
			{
				Wlist <- InitialMatrix(X, y, ycat, nslices);

				values <- c(objfun(W=list(Qt=Wlist$Wsir, dim=c(d,p), Sigmas=Sigmas))$value, 
					objfun(W=list(Qt=Wlist$Wsave, dim=c(d,p), Sigmas=Sigmas))$value, 
					objfun(W=list(Qt=Wlist$Wdr, dim=c(d,p), Sigmas=Sigmas))$value);

				W <- list(Qt=Wlist[[which(values == max(values))[1]]], dim=c(d,p), Sigmas=Sigmas)
			}
			if (!is.null(mf$sim_anneal)) W <- list(dim=c(d,p), Sigmas=Sigmas);

			grassmann <- GrassmannOptim(objfun, W, verbose=verbose,...);

			Gammahat = matrix(grassmann$Qt[,1:d], ncol=d);

			invDeltahat <- solve(Sigmatilde) + 
				Gammahat%*%solve(t(Gammahat)%*%Deltatilde%*%Gammahat)%*%t(Gammahat) -
				Gammahat%*%solve(t(Gammahat)%*%Sigmatilde%*%Gammahat)%*%t(Gammahat);

			Deltahat <- solve(invDeltahat); Pj <- projection(Gammahat, Deltahat);

			Deltahat_y <- vector("list", nslices);

			for (i in 1:nslices)
			{
				Deltahat_y[[i]] <- Deltahat + t(Pj) %*% (Deltatilde_y[[i]] - Deltahat)%*%Pj;
			}
			terme0 <- -(n*p)*(1+log(2*pi))/2;

			terme1 <- -n*log( det(Sigmatilde) )/2; 

			terme2 <- n*log(det(t(Gammahat)%*%Sigmatilde%*%Gammahat))/2;

			terme3<-0;

			for (i in 1:length(N_y))
			{
				terme3 <- terme3 - N_y[i] * log( det(t(Gammahat)%*%Deltatilde_y[[i]]%*%Gammahat) )/2;
			} 
			loglik = terme0 + terme1 + terme2 + terme3;

			numpar <- p+(nslices-1)*d+p*(p+1)/2+d*(p-d)+(nslices-1)*d*(d+1)/2; 

			AIC <- -2*loglik + 2 * numpar;

			BIC <- -2*loglik + log(n) * numpar;
		}
		return(list(Gammahat=Gammahat, Deltahat=Deltahat, Deltahat_y=Deltahat_y, loglik=loglik, 
				numpar=numpar, aic=AIC, bic=BIC))
	}

	if (identical(short, TRUE))
	{
		fit <- onelad(d);

		ans <- list(Gammahat=fit$Gammahat, Deltahat=fit$Deltahat, Deltahat_y=fit$Deltahat_y, loglik=fit$loglik, 
			aic=fit$aic, bic=fit$bic, numpar=fit$numpar, numdir=d, model="lad", call=match.call(expand.dots=TRUE),short=short);

		class(ans) <- "ldr";	
		
		return(invisible(ans));	
	}
	aic <- bic <- numpar <- loglik <- vector(length=d+1);

	Deltahat <- Deltahat_y <- Gammahat <-vector("list");  
			
	for (i in 0:d)
	{
		fit <- onelad(i);
			
		Gammahat[[i+1]] <-fit$Gammahat; 

		Deltahat[[i+1]] <- fit$Deltahat;

		Deltahat_y[[i+1]] <- fit$Deltahat_y;

		loglik[i+1] <- fit$loglik; 

		numpar[i+1] <- fit$numpar;

		aic[i+1] <- fit$aic; 

		bic[i+1] <- fit$bic;	
	}
	ans <- list(Gammahat=Gammahat, Deltahat=Deltahat, Deltahat_y=Deltahat_y, loglik=loglik, aic=aic, 
			bic=bic, numpar=numpar, numdir=d, model="lad", call=match.call(expand.dots=TRUE), short=short);

	class(ans)<- "ldr";

	return(invisible(ans))
}
