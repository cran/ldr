pfc <-
function(X, fy=fy, numdir=2, structure=c("iso", "aniso", "unstr", "unstr2"), 
    eps_aniso=1e-2, verbose=FALSE, short=TRUE,...)
{
	if (numdir <= 0){cat("`numdir` does not have a valid entry"); return()}

	"%^%"<-function(M, p) 
	{ 
		# This operator compute the M^p where M is a matrix;

		if (!is.matrix(M)) stop("Argument of power is not a matrix in pfc")

		if (prod(dim(M)==list(1,1))) return( as.matrix(M^p) );

		svdM <- svd(M); return(svdM$u%*%diag(c(svdM$d)^p)%*%t(svdM$v));
	}

	Trace<-function(X)
	{	
		if (!is.matrix(X)) stop("Argument to Trace is not a matrix in pfc");

		return(sum(diag(X))) 
	}
	onepfc = function(X, fy)
	{
		# X is univariate predictor

		nobs <- length(X); r <- dim(fy)[2]; 

		P_F <- fy%*%solve(t(fy)%*%fy)%*%t(fy);

		Xc <- scale(X, TRUE, FALSE)

		Sigmahat_fit <- (1/nobs)*t(Xc)%*%P_F%*%(Xc);

		ev.fit <- eigen(Sigmahat_fit);


		temp.dat<-data.frame(cbind(X, fy)); xnam<-paste("xx", 1:r, sep="");

		names(temp.dat)<-c("yy", xnam);

		fm.lm<- as.formula( paste("yy ~ ", paste(xnam, collapse= "+")));

		summary.fm <- summary(lm(fm.lm, data=temp.dat));

		Betahat <- matrix(summary.fm$coefficients[2:(r+1),1], ncol=r);

		Gammahat <- matrix(1, ncol=1, nrow=1); 

		Deltahat <- matrix(summary.fm$sigma^2, ncol=1, nrow=1);

		Muhat <- matrix(summary.fm$coefficients[1,1], ncol=1);


		loglik <- - 0.5*n*(1+log(2*pi*summary.fm$sigma^2));

		numpar <- p  + dim(fy)[2] + 1;

		aic <- -2*loglik + 2*numpar;

		bic <- -2*loglik + log(n)*numpar;

		ans <- list(Muhat=Muhat, Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, 
				loglik=loglik, aic=aic, bic=bic, numpar=numpar, numdir=1, model="pfc", 
				call=match.call(), structure="iso", fy=fy, Xc=Xc, short=short);

		class(ans)<- "ldr";

		return(ans)
	}
	r <- dim(fy)[2]; X <- as.matrix(X); 

	op <- dim(X); n <- op[1]; p <- op[2]; 

	eff.numdir <- min(numdir, r, p);	

	vnames <- dimnames(X)[[2]];

	if (is.null(vnames)) vnames <- paste("X", 1:p, sep="");

	if (p==1) return(onepfc(X, fy));
	

	Muhat <- apply(X, 2, mean);

	Xc <-  scale(X, TRUE, FALSE);

	P_F <- fy%*%solve(t(fy)%*%fy)%*%t(fy);

	Sigmahat <- t(Xc)%*%Xc/n;

	Sigmahat_fit <- t(Xc)%*%P_F%*%Xc/n; 

	Sigmahat_res <- Sigmahat - Sigmahat_fit;


	if (structure=="iso")
	{
		iso <- function(i)
		{
			ev <- eigen(Sigmahat);	

			ev.fit <- eigen(Sigmahat_fit); 

			all_evalues <-ev.fit$values

			evalues <- all_evalues[1:i]

			sigma2hat <- Re(sum(ev$values)/p); 

			Gammahat <- Re(matrix(ev.fit$vectors[,1:i], ncol=i));

			dimnames(Gammahat) <- list(vnames, paste("Dir", 1:i, sep=""))

			Betahat <-Re(t(Gammahat)%*%t(Xc)%*%fy%*%solve(t(fy)%*%fy));

			sigma2hat <- Re((sum(ev$values)-sum(evalues))/p); 

			Deltahat <- sigma2hat*diag(1, p); 

			dimnames(Deltahat) <- list(vnames, vnames)

			loglik <- - 0.5*n*p*(1+log(2*pi*sigma2hat));

			numpar <- p + (p-i)*i + i*dim(fy)[2] + 1;

			aic <- -2*loglik + 2*numpar;

			bic <- -2*loglik + log(n)*numpar;

			return(list(Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, 
					evalues=evalues, loglik=loglik, aic=aic, bic=bic, numpar=numpar));
		}

		if (identical(short, TRUE))
		{
			out <- iso(eff.numdir);

			ans <- list(Muhat=Muhat, Betahat=out$Betahat, Gammahat=out$Gammahat, 
					Deltahat=out$Deltahat, loglik=out$loglik, aic=out$aic, 
					bic=out$bic, numpar=out$numpar, numdir=eff.numdir, model="pfc", 
					evalues=out$evalues, structure="iso", fy=fy,  Xc=Xc, 
					call=match.call(expand.dots=TRUE), short=short);

			class(ans) <- "ldr";	
		
			return(ans);		
		}

		if (identical(short, FALSE))
		{
			aic <- bic <- numpar <- loglik <- vector(length=eff.numdir+1);

			Betahat <- Deltahat <- Gammahat <-vector("list"); 

			# No fitting values (eff.numdir=0)

			ev <- eigen(Sigmahat); 

			sigma2hat <- sum(ev$values)/p; 

			loglik[1] <- - 0.5*n*p*(1+log(2*pi*sigma2hat));

			numpar[1] <- p + 1;

			aic[1] <- -2*loglik[1] + 2*numpar[1];

			bic[1] <- -2*loglik[1] + log(n)*numpar[1];
		

			for (i in 1:eff.numdir)
			{
				fit <- iso(i);
			
				Betahat[[i]] <-fit$Betahat; 

				Gammahat[[i]] <-fit$Gammahat; 

				Deltahat[[i]] <- fit$Deltahat;

				loglik[i+1] <- fit$loglik; 

				numpar[i+1] <- fit$numpar;

				aic[i+1] <- fit$aic; 

				bic[i+1] <- fit$bic;	
			}
			ans <- list(Muhat=Muhat, Betahat=Betahat, Gammahat=Gammahat, 
					Deltahat=Deltahat, loglik=loglik, aic=aic, bic=bic, numpar=numpar, 
					numdir=eff.numdir, model="pfc", evalues=fit$evalues, structure="iso", 
					fy=fy,  Xc=Xc, call=match.call(), short=short);

			class(ans)<- "ldr";

			return(ans)
		} 
	}

	if (structure=="aniso")
	{
		aniso = function(X, d, fy, eps_aniso=1e-4, short)
		{
			vnames <- dimnames(X)[[2]];

			if (is.null(vnames)) vnames <- paste("X", 1:ncol(X), sep="");

			op <- dim(X); n <- op[1]; p <- op[2];

			# Initial Step
			fit <- pfc(X=X, numdir=d, fy=fy, structure="iso", short=short);

			if (identical(short, TRUE))
			{
				Betahatx <- fit$Betahat; Gammahatx <- fit$Gammahat; 

				Xc <- scale(X, TRUE, FALSE) - fy%*%t(Gammahatx%*%Betahatx);

				deltahat <- diag(t(Xc)%*%(Xc)/n);

				repeat
				{
					Xnew = X%*%((1/sqrt(deltahat))*diag(p));

					fit <- pfc(X=Xnew, numdir=d, fy=fy, structure="iso", short=TRUE);

					Betahatx <- fit$Betahat; Gammahatx <- (diag(p)*sqrt(deltahat))%*%fit$Gammahat; 

					Xc <- scale(X, TRUE, FALSE) - fy%*%t(Gammahatx%*%Betahatx);

					deltahat0 <- diag(t(Xc)%*%(Xc)/n);

					if (sum(abs(deltahat-deltahat0)) < (eps_aniso*p)) break;

					deltahat <- deltahat0; 
				}
				dimnames(Gammahatx) <- list(vnames, paste("Dir", 1:d, sep=""))

				Deltahat <- deltahat*diag(p);

				dimnames(Deltahat) <- list(vnames, vnames);

				loglik <- - 0.5*n*p*(1+log(2*pi)) - 0.5*n*log(prod(deltahat));

				numpar <- p + d*(p-d) + ncol(fy)*d + p; 

				aic <- -2*loglik + 2*numpar;

				bic <- -2*loglik + log(n)*numpar;

				ans <- list(Betahat=Betahatx, Gammahat=OrthoNorm(Gammahatx), Deltahat=Deltahat, 
						evalues=fit$evalues, loglik=loglik, aic=aic, bic=bic, numpar=numpar, short=short);
		
				return(ans)
			}

			Deltahat <- Betahat <- Gammahat <- vector("list");

			aic <- bic <- numpar <- loglik <- vector(length=eff.numdir+1);

			# No fitting values (eff.numdir=0)

			Xc <- scale(X, TRUE, FALSE);

			Sigmahat <- t(Xc)%*%Xc/n;

			ev <- eigen(Sigmahat); 

			loglik[1] <- - 0.5*n*p*(1+log(2*pi)) - 0.5*n*log(prod(diag(Sigmahat)));

			numpar[1] <- p + p;

			aic[1] <- -2*loglik[1] + 2*numpar[1];

			bic[1] <- -2*loglik[1] + log(n)*numpar[1];


			for (i in 1:eff.numdir)
			{
				Betahatx <- fit$Betahat[[i]]; Gammahatx <- fit$Gammahat[[i]]; 

				Xc <- scale(X, TRUE, FALSE) - fy%*%t(Gammahatx%*%Betahatx);

				deltahat <- diag(t(Xc)%*%(Xc)/n);

				repeat
				{
					Xnew = X%*%((1/sqrt(deltahat))*diag(p));

					fit2 <- pfc(X=Xnew, numdir=i, fy=fy, structure="iso", short=TRUE);

					Betahatx <- fit2$Betahat; Gammahatx <- (diag(p)*sqrt(deltahat))%*%fit2$Gammahat; 

					Xc <- scale(X, TRUE, FALSE) - fy%*%t(Gammahatx%*%Betahatx);

					deltahat0 <- diag(t(Xc)%*%(Xc)/n);

					if (sum(abs(deltahat-deltahat0)) < (eps_aniso*p)) break;

					deltahat <- deltahat0; 
				}

				Deltahat[[i]] <- deltahat*diag(p);

				dimnames(Deltahat[[i]]) <- list(vnames, vnames);

				temp0 <- -(n*p/2)*(1 + log(2*pi));

				temp1 <- -(n/2)*sum( log( eigen(Deltahat[[i]], symmetric=T)$values) ); 

				temp2 <- 0;

				if (i < min(ncol(fy), p))
				{
					Xc <-  scale(X, TRUE, FALSE);	P_F <- Re(fy%*%solve(t(fy)%*%fy)%*%t(fy));

					Sigmahat <- t(Xc)%*%Xc/n; Sigmahat_fit <- t(Xc)%*%P_F%*%Xc/n; 

					Sigmahat_res <- Sigmahat - Sigmahat_fit;

					sqrt_Sigmahat_res <- Re(Sigmahat_res%^%0.5); 

					Inv_Sqrt_Sigmahat_res <- solve(sqrt_Sigmahat_res);

					lf_matrix <- Inv_Sqrt_Sigmahat_res%*% Sigmahat_fit %*% Inv_Sqrt_Sigmahat_res;

					temp2 <- -(n/2)*sum(log(1+ eigen(lf_matrix, symmetric=T)$values[(i+1):min(ncol(fy), p)]));
				}
				loglik[i+1] = temp0 + temp1 + temp2;

				numpar[i+1] <- p + (p-i)*i + i*dim(fy)[2] + p;

				aic[i+1] <- -2*loglik[i+1] + 2*numpar[i+1];

				bic[i+1] <- -2*loglik[i+1] + log(n)*numpar[i+1];

				Betahat[[i]] <- Betahatx;

				Gammahat[[i]] <- OrthoNorm(Gammahatx);

				dimnames(Gammahat[[i]]) <- list(vnames, paste("Dir", 1:i, sep=""))
			}
			ans <- list(Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, evalues=fit2$evalues, 
						loglik=loglik, aic=aic, bic=bic, numpar=numpar, short=short);

			return(ans)	
		}

		fit <- aniso(X=X, d=eff.numdir, fy=fy, eps_aniso=eps_aniso, short=short);

		ans <- list(Muhat=Muhat, Betahat=fit$Betahat, Gammahat=fit$Gammahat, Deltahat=fit$Deltahat, model="pfc",
				loglik=fit$loglik, aic=fit$aic, bic=fit$bic, numpar=fit$numpar, numdir=eff.numdir, 
				evalues=fit$evalues, structure="aniso", Xc=Xc, fy=fy,  call=match.call(), short=fit$short);

		class(ans)<- "ldr";

		return(ans)
	}


	if (structure=="unstr")
	{
		unstr<-function(i)
		{
			sqrt_Sigmahat_res <- Sigmahat_res%^%0.5; 

			Inv_Sqrt_Sigmahat_res <- solve(sqrt_Sigmahat_res);

			lf_matrix <- Inv_Sqrt_Sigmahat_res%*%Sigmahat_fit%*%Inv_Sqrt_Sigmahat_res;

			all_evalues <- eigen(lf_matrix, symmetric=T)$values;

			evalues <- all_evalues[1:i];


			Vhat <- eigen(lf_matrix, symmetric=T)$vectors;

			Vhati <- matrix(Vhat[,1:i], ncol=i);

			Gammahat <- (Sigmahat_res%^%0.5)%*%Vhati%*%solve((t(Vhati)%*%Sigmahat_res%*%Vhati)%^%0.5);  

			dimnames(Gammahat)<- list(vnames, paste("Dir", 1:i, sep=""));


			Khat<-diag(0, p); 

			if (i < min(ncol(fy),p)) {diag(Khat)[(i+1):min(ncol(fy), p )]<- all_evalues[(i+1):min(ncol(fy), p)]};

			Deltahat <- sqrt_Sigmahat_res%*%Vhat%*%(diag(p)+Khat)%*%t(Vhat)%*%sqrt_Sigmahat_res;

			dimnames(Deltahat) <- list(vnames, vnames);

			Betahat <- ((t(Vhati)%*%Sigmahat_res%*%Vhati)%^%0.5)%*%t(Vhati)%*%solve(Sigmahat_res%^%0.5)%*%t(Xc)%*%fy%*% solve(t(fy)%*%fy);

			temp0 <- -(n*p/2)*(1 + log(2*pi));

			temp1 <- -(n/2)*sum(log(svd(Sigmahat_res)$d[1:i])); 

			temp2 <- 0; 

			if (i < min(ncol(fy),p)) temp2 <- -(n/2)*sum(log(1+ all_evalues[(i+1):min(ncol(fy), p)]));

			loglik <- temp0 + temp1 + temp2;

			numpar <- p + (p-i)*i + i*ncol(fy) + p*(p+1)/2;

			aic <- -2*loglik + 2*numpar;

			bic <- -2*loglik + log(n)*numpar;

			return(list(Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, evalues=evalues, loglik=loglik, aic=aic, bic=bic, numpar=numpar));
		}

		if (identical(short, TRUE))
		{
			out <- unstr(eff.numdir);

			ans <- list(Muhat=Muhat, Betahat=out$Betahat, Gammahat=out$Gammahat, Deltahat=out$Deltahat, 
					evalues=out$evalues, loglik=out$loglik, aic=out$aic, bic=out$bic, numpar=out$numpar, 
					numdir=eff.numdir, model="pfc", structure="unstr", fy=fy, Xc=Xc,  call=match.call(), short=short);

			class(ans) <- "ldr";	
		
			return(ans);	
		}

		aic <- bic <- numpar <- loglik <- vector(length=eff.numdir+1);

		evalues <- vector(length=eff.numdir);

		Betahat <- Deltahat <- Gammahat <-vector("list"); 
 
		loglik[1] <- - 0.5*n*p*(1+log(2*pi)) - 0.5*n*log(det(Sigmahat)) ;

		numpar[1] <- p + p*(p+1)/2; aic[1] <- -2*loglik[1] + 2*numpar[1];

		bic[1] <- -2*loglik[1] + log(n)*numpar[1];

		Deltahat[[1]] <- Sigmahat; dimnames(Deltahat[[1]]) <- list(vnames, vnames); 

		for (i in 1:eff.numdir)
		{
			fit <- unstr(i);
			
			Betahat[[i]] <-fit$Betahat; 

			Gammahat[[i]] <-fit$Gammahat; 

			Deltahat[[i]] <- fit$Deltahat;

			loglik[i+1] <- fit$loglik; 

			numpar[i+1] <- fit$numpar;

			aic[i+1] <- fit$aic; 

			bic[i+1] <- fit$bic;	
		}
		ans <- list(Muhat=Muhat, Betahat=Betahat, Gammahat=Gammahat, Deltahat=Deltahat, 
				evalues=fit$evalues, loglik=loglik, aic=aic, bic=bic, numpar=numpar, 
				numdir=eff.numdir, model="pfc", structure="unstr", fy=fy, Xc=Xc,  call=match.call(), short=short);

		class(ans)<- "ldr";

		return(ans)
	}

	else if (structure=="unstr2")
	{
		unstr2 <- function(i)
		{
			objfun <- function(W)
			{
				Qt <- W$Qt; dc <- W$dim[1]; 

				p <- ncol(Qt); S <- W$Sigmas;
 
				U <- matrix(Qt[,1:dc], ncol=dc);	

				V <- matrix(Qt[,(dc+1):p], ncol=(p-dc));

				value <- -(n/2)*(p*log(2*pi)+p+log(det(t(V)%*%S$Sigmahat%*%V))+ log(det(t(U)%*%S$Sigmahat_res%*%U)));

				terme1 <- solve(t(U)%*%S$Sigmahat_res%*%U)%*%(t(U)%*%S$Sigmahat_res%*%V); 

				terme2 <- (t(U)%*%S$Sigmahat%*%V)%*%solve(t(V)%*%S$Sigmahat%*%V);

				gradient <- 2*(terme1 - terme2);

				return(list(value=value, gradient=gradient));
			}
			sigmas <- list(Sigmahat=Sigmahat, Sigmahat_fit=Sigmahat_fit, Sigmahat_res=Sigmahat_res, p=p, n=n);

			W <- list(Qt = svd(Sigmahat_fit)$u, dim=c(numdir, p), Sigmas=list(Sigmahat=Sigmahat, Sigmahat_fit=Sigmahat_fit, 
						Sigmahat_res=Sigmahat_res, p=p, n=n));

			objfun <- assign("objfun", objfun, envir=.GlobalEnv); 

			grassoptim <- GrassmannOptim(objfun, W, verbose=verbose, ...);

			Gammahat <- matrix(grassoptim$Qt[,1:i], ncol=i, dimnames=list(vnames, 
						paste("Dir", 1:i, sep=""))); 

			Gammahat0 <- matrix(grassoptim$Qt[, (i+1):p], ncol=p-i, dimnames=list(vnames, paste("Dir", 
						(i+1):p, sep="")));

			Betahat <- t(Gammahat)%*%t(Xc)%*%fy%*%solve(t(fy)%*%fy);

			Omegahat <- t(Gammahat)%*%Sigmahat_res%*%Gammahat;

			Omegahat0 <-t(Gammahat0)%*%Sigmahat%*%Gammahat0;

			Deltahat <- Gammahat%*%Omegahat%*%t(Gammahat) + Gammahat0%*%Omegahat0%*%t(Gammahat0);

			dimnames(Deltahat) <- list(vnames, vnames);


			temp0 <- -(n*p/2)*(1+log(2*pi));

			temp1 <- -(n/2)*log(det(t(Gammahat)%*%Sigmahat_res%*%Gammahat)); 

			temp2 <- -(n/2)*log(det(t(Gammahat0)%*%Sigmahat%*%Gammahat0)); 

			loglik <- temp0 + temp1 + temp2;  

			numpar <- p + (p-i)*i + i*dim(fy)[2] + i*(i+1)/2 + (p-i)*(p-i+1)/2; 

			aic <- -2*loglik + 2*numpar;

			bic <- -2*loglik + log(n)*numpar;

			ev.fit <- eigen(Sigmahat_fit); 

			evalues <- ev.fit$values[1:i];

			return(list(Betahat=Betahat, Gammahat=Gammahat, Gammahat0=Gammahat0, Omegahat=Omegahat, Omegahat0=Omegahat0, 
					Deltahat=Deltahat, evalues=evalues, loglik=loglik, aic=aic, bic=bic, numpar=numpar));
		}

		if (identical(short, TRUE))
		{
			out <- unstr2(numdir);

			ans <- list(Muhat=Muhat, Betahat=out$Betahat, Gammahat=out$Gammahat, Gammahat0=out$Gammahat0, 
					Omegahat=out$Omegahat, Omegahat0=out$Omegahat0, Deltahat=out$Deltahat, 
					evalues=out$evalues, loglik=out$loglik, aic=out$aic, bic=out$bic, 
					numpar=out$numpar, numdir=numdir, model="pfc", structure="unstr2", 
					fy=fy, Xc=Xc,  call=match.call(), short=short);

			class(ans) <- "ldr";	
		
			return(ans);
		}

		aic <- bic <- numpar <- loglik <- vector(length=numdir+1);

		Betahat <- Deltahat <- Gammahat <- Gammahat0 <- Omegahat <- Omegahat0 <- vector("list"); 

		loglik[1] <- -(n*p/2)*(log(2*pi) + (1+log(Trace(Sigmahat)/p)));

		numpar[1] <- p + p*(p+1)/2;

		aic[1] <- -2*loglik[1] + 2*numpar[1];

		bic[1] <- -2*loglik[1] + log(n)*numpar[1];

		for(m in 1:numdir)
		{
			fit <- unstr2(m);

			Betahat[[m]] <-fit$Betahat; 

			Gammahat[[m]] <-fit$Gammahat; 
	
			Omegahat[[m]] <- fit$Omegahat; 

			Omegahat0[[m]] <- fit$Omegahat0;

			Deltahat[[m]] <- fit$Deltahat; 

			loglik[m+1] <- fit$loglik; 

			numpar[m+1] <- fit$numpar;

			aic[m+1] <- fit$aic; 

			bic[m+1] <- fit$bic;	
		}
		ans <- list(evalues=fit$evalues, loglik =loglik, aic=aic, bic=bic, numdir=numdir, 
				numpar=numpar, Muhat=Muhat, Betahat=Betahat, Gammahat=Gammahat,	
				Gammahat0=Gammahat0, Omegahat=Omegahat, Omegahat0=Omegahat0, 
				Deltahat=Deltahat, model="pfc", structure="unstr2", 
				fy=fy, Xc=Xc,  call=match.call(), short=short);

		class(ans)<- "ldr";

		return(ans)
	}

}
