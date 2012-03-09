library(geepack)
library(Matrix)
library(lme4)
library(ncdf)

updateR <- function (xnew, R = NULL, xold, eps = .Machine$double.eps){
    R <- as.matrix(R)
    xtx <- as.matrix(crossprod(xnew))
    norm.xnew <- sqrt(xtx)
    if (is.null(R)) {
        R <- matrix(norm.xnew, 1, 1)
        attr(R, "rank") <- 1
        return(R)
    }
    Xtx <- as.numeric(t(xnew) %*% xold)
    k <- NCOL(R)
    nb <- as.integer(1)
    job <- as.integer(11)
    r <- .C("bakslv", t = R, ldt = k, n = k, b = Xtx, ldb = k, 
        nb = nb, x = matrix(0, k, nb), job = job, info = integer(1L), 
        DUP = FALSE, PACKAGE = "base")["x"]$x
    rpp <- xtx - crossprod(r)
    rank <- attr(R, "rank")
    rpp <- sqrt(rpp)
    rank <- rank + 1
    R <- cbind(rbind(as.matrix(R), 0), c(r, rpp))
    attr(R, "rank") <- rank
    R
}

updateR2 <- function (xnew, R = NULL, xold, eps = .Machine$double.eps){
    R <- as.matrix(R)
    xtx <- crossprod(xnew)
    if (is.null(R)) {
    	norm.xnew <- chol(xtx)
        R <- matrix(norm.xnew, 1, 1)
        attr(R, "rank") <- 1
        return(R)
    }
    Xtx <- as.matrix(t(xnew) %*% xold)
    k <- NCOL(R)
    nb <- as.integer(NCOL(xnew))
    job <- as.integer(11)
    r <- .C("bakslv", t = R, ldt = k, n = k, b = t(Xtx), ldb = k, 
        nb = nb, x = matrix(0, k, nb), job = job, info = integer(1L), 
        DUP = FALSE, PACKAGE = "base")["x"]$x
    rpp <- chol(xtx - crossprod(r))
    rank <- attr(R, "rank")
    rank <- rank + 2
    R <- cbind(rbind(R, 0,0), rbind(r, rpp))
    attr(R, "rank") <- rank
    R
}

genCor <- function(maxClustSize, alpha, corstr){
	if(corstr=="independence"){#indep
		return(as.matrix(Diagonal(x=rep(1,maxClustSize))))	
	}
	if(corstr=="exchangeable"){#exch
			return(alpha^(!outer(1:maxClustSize,1:maxClustSize,"==")))
	}
	if(corstr=="ar1"){#ar1
			return(alpha^abs(outer(1:maxClustSize,1:maxClustSize,"-")))
	}
	stop("Unrecognized correlation structure, should be one of 'independence', 'exchangeable' or 'ar1'")
}

boss.set <- function(formula, E.name=NULL, family=gaussian(), id = NULL, corstr = "independence", type="glm", method="chol", data,...){#
	init <- NULL
									
	init$formula <- formula						
	init$family <- family
	init$method = method
	m <- glm(formula,family=family, data=data)
	if(is.null(id)){
		id <- 1:NROW(m$res)
	}
	init$id <- id
	##subset to complete cases
	if(!is.null(m$na.action)){
		data <- data[-m$na.action,]
		init$id <- id[-m$na.action]
	}

	init$y <- m$model[,1]
	init$X <- model.matrix(m)
	init$p <- NCOL(init$X)
	init$n <- NROW(init$X)
	
	if(type == "glm"){ ##independet observations
		init$w <- family$variance(m$fit)
		init$d <- sqrt(init$w)
		init$yhat <- m$fit
		Wdi <- Diagonal(init$n,sqrt(init$w))
		init$b <- coef(m)
		init$XD <- Wdi%*%init$X
		init$R <- chol(crossprod(init$XD))
		
		if(!is.null(E.name)){
			if( !(E.name %in% attr(terms(formula),"term.labels")) ){
				stop("Main effect of E must be included")
			}
			init$E <- data[,E.name]
			init$Ed <- Wdi%*%init$E
			init$Ev <- init$E
		}
		if(family$family=="gaussian"){
			init$XY <- crossprod(init$X,init$y)
			if(!is.null(E.name) ) class(init) <- "flm.GxE"			else class(init) <- "flm"
		} else if(family$family=="binomial"){
			init$XY <- crossprod(init$X,init$y-init$yhat)
			if(!is.null(E.name) ) class(init) <- "flr.GxE"			else class(init) <- "flr"	
		} else {
			stop("family not supported")	
		}			
	} else if(type == "gee"){ ##GEE
		m <- geese(formula, family=family, id= id, corstr=corstr, data = data)
		clust.sizes <- m$clusz
	
		init$corstr <- corstr
		init$yhat <- m$fit <- family$linkinv(init$X%*%m$beta)
		init$w <- family$variance(m$fit)
		init$b <- coef(m)
			
		V <- 1
		Vi <- 1
			
		m$working.correlation <- genCor(max(clust.sizes),m$alpha,corstr)
		if(max(clust.sizes) > 1){	
			for(i in 2:max(clust.sizes)){
				V <- c(V,list(m$working.correlation[1:i,1:i]))
				Vi <- c(Vi,list(solve(V[[ i ]])))
			}
		}
		Vi.sqrt <- lapply(Vi, chol) 
	
		Wd <- Diagonal(init$n, 1/sqrt(init$w))
		Wdi <- Diagonal(init$n,sqrt(init$w))
		big.V <- bdiag(lapply(clust.sizes,function(k){V[[k]]}))
		big.Vi <- bdiag(lapply(clust.sizes,function(k){Vi[[k]]}))
		big.Vi.sqrt <- bdiag(lapply(clust.sizes,function(k){Vi.sqrt[[k]]}))
		init$XD <- as.matrix((big.Vi.sqrt%*%Wdi)%*%init$X)
		init$XVi <- as.matrix((Wd%*%big.Vi%*%Wdi)%*%init$X)
		init$R <- as.matrix(chol(crossprod(init$XD)))
	
		init$v <- rowSums(Wd%*%big.Vi%*%Wd)
		init$d <- rowSums(Wd%*%big.Vi.sqrt)
		
		init$vw <- rowSums(Wd%*%big.Vi%*%Wdi)
		init$dw <- rowSums(big.Vi.sqrt%*%Wdi)
	
		if(!is.null(E.name)){
			if( !(E.name %in% attr(terms(formula),"term.labels")) ){
				stop("Main effect of E must be included")
			}
			init$E <- init$Ev <- init$Ed <- data[,E.name]
			init$Ed <- as.numeric(big.Vi.sqrt%*%Wdi%*%init$E)
			init$Ev <- as.numeric(Wd%*%big.Vi%*%Wdi%*%init$E)
		}
		if(family$family=="gaussian"){
			init$XY <- crossprod(init$XVi,init$y)
			if(!is.null(E.name) ) class(init) <- "fgee.GxE"			else class(init) <- "fgee"
		} else if(family$family=="binomial"){
			init$XY <- crossprod(init$XVi,init$y-init$yhat)
			if(!is.null(E.name) ) class(init) <- "fgee.lr.GxE"			else class(init) <- "fgee.lr"	
		} else {
			stop("family not supported")	
		}
		
	} else if(type =="lmm"){ #linear mixed models
		if(init$family$family != "gaussian") stop("Only linear mixed models supported at this time")
			
		m <- lmer(formula, data = data)
		init$X <- m@X
		init$p <- ncol(init$X)
		init$yhat <- m@mu
		init$data <- data
		
		V <- crossprod(m@A)+Diagonal(init$n,1)
		V.sqrt <- chol(V)
		Vi.sqrt <- solve(V.sqrt)
		Vi <- Vi.sqrt%*%t(Vi.sqrt)
		init$Vi.sqrt <- Vi.sqrt
		init$Vi <- Vi
		init$XD <- as.matrix(t(Vi.sqrt)%*%init$X)
		init$XVi <- as.matrix(Vi%*%init$X)
		init$R <- chol(crossprod(init$XD))
		init$XY <- t(init$XVi)%*%init$y
		
		if(!is.null(E.name)){
			if( !(E.name %in% attr(terms(formula),"term.labels")) ){
				stop("Main effect of E must be included")
			}
			init$E <- init$Ev <- init$Ed <- data[,E.name]
			init$Ed <- as.numeric(t(Vi.sqrt)%*%init$E)
			init$Ev <- as.numeric(Vi%*%init$E)
			
			class(init) <- "lmm.GxE"
		} else {
			class(init) <- "lmm"
		}
	}	
	
	#Approximations for single marker models:	
	if(method=="swap"){
		if(!is.null(E.name)) stop("Only single marker models supported with 'swap' method")
		
		if( type == "lmm"){#Linear Mixed models
			
			init$Wdi <- t(Vi.sqrt)
			init$z <- as.numeric(init$Wdi%*%init$y)
			
			if( all(  round(t(Vi.sqrt)%*%as.numeric(init$id) - rowSums(t(Vi.sqrt))*as.numeric(init$id),6)==0) ){
				#genotype is constant in independent clusters
				init$d <- rowSums(init$Wdi)
				class(init) <- "smcg"
			} else {
				#genotype varies in clusters
				class(init) <- "smvg"	
			} 
				
		} else if(type=="glm"){  #independent observations
			
			init$z <- init$d*(init$X%*%init$b + (init$y-init$yhat)/init$w)
			class(init) <- "smcg"
			
		} else if(type=="gee"){ #GEE with non-robust standard errors:
			
			init$Wdi <- big.Vi.sqrt%*%Wdi
			init$d <- init$dw
			#init$z <- as.numeric(init$Wdi%*%sqrt(init$w)*(init$yhat + (init$y-init$yhat)/init$w))
			init$z <- as.numeric(init$Wdi%*%(init$yhat + (init$y-init$yhat)/init$w))
			class(init) <- "smcg"
		}
		init$Xy <- cbind(as.matrix(init$XD),init$z)
		init$A <- solve(crossprod(init$Xy))		
		init$AXy <- init$A%*%t(init$Xy)
		init$yres <- var(as.numeric(init$z - init$XD%*%(solve(crossprod(init$XD))%*%t(init$XD)%*%init$z)))
	}	
	
	return(init)
}

boss.fit <- function(g, init, thresh = 1e-7, robust=TRUE,...){ 
	UseMethod("boss.fit", init)
}

boss.fit.smcg <- function(g,init, thresh = 1e-7, robust= TRUE,...){
	g <- as.matrix(g)
	dg <- as.matrix(init$d*g)
	betas <- init$AXy%*%dg
	res <- dg - init$Xy%*%betas 
	res.vars <- colSums(res^2)/(init$n-init$p-1)
	chi2s <- betas[init$p+1,]^2/(res.vars*init$A[init$p+1,init$p+1])
	
	#set missing those regressions where colinearity would prevent matrix inversion in standard methods:
	chi2s[zapsmall(res.vars) == 0] <- NA	
	
	beta.main <- betas[init$p+1, ]*(init$yres/res.vars)
	v.main <- beta.main^2/chi2s
	
	return(list(beta.main = beta.main,v.main = v.main, chi2 = chi2s))
	}
	
boss.fit.smvg <- function(g,init,thresh = 1e-7, robust= TRUE,...){
	g <- as.matrix(g)
	dg <- init$Wdi*g
	betas <- init$AXy%*%dg
	res <- dg - init$Xy%*%betas 
	res.vars <- colSums(res^2)/(init$n-init$p)
	chi2s <- betas[init$p,]^2/(res.vars*init$A[init$p,init$p])
	
	#set missing those regressions where colinearity would prevent matrix inversion in standard methods:
	chi2s[zapsmall(res.vars) == 0] <- NA	
	
	beta.main <- betas[init$p, ]*(init$yres/res.vars)
	v.main <- chi2s/beta.main^2
	
	return(list(beta.main = beta.main, v.main = v.main, chi2s = chi2s))
	}
	
boss.fit.lmm <- function(g,init, thresh = 1e-7, robust = TRUE, ...){
	g <- as.numeric(g)
	p <- init$p+1
	
	R <- updateR(as.matrix(t(init$Vi.sqrt)%*%g),init$R,init$XD)
	X <- cbind(init$X,g)
	
	A <- chol2inv(R)
	beta <- A%*%c(init$XY,sum(as.matrix(init$Vi%*%g)*init$y))
	
	r <- init$y-X%*%beta
	V <- sum(r^2)/(init$n-p)*A	
	
	if(1-pchisq(beta[p]^2/V[p,p],1) < thresh){
		mod <- lmer(update(init$formula,.~.+g), data = init$data,...)
		beta <- mod@fixef
		V <- vcov(mod)
		k <- match("g",names(beta))
		return(list(beta.main = beta[k], v.main = V[k,k]))		}
	
	return(list(beta.main=beta[p], v.main = V[p,p]))
	}
	
boss.fit.lmm.GxE <- function(g,init, thresh = 1e-7, robust= TRUE,...){
	g <- as.numeric(g)
	p <- init$p+2
	
	R <- updateR2(cbind(as.numeric(t(init$Vi.sqrt)%*%g),g*init$Ed),init$R,init$XD)
	X <- cbind(init$X,g,g*init$E)
	
	A <- chol2inv(R)
	beta <- A%*%c(init$XY,sum(as.matrix(init$Vi%*%g)*init$y),sum(as.matrix(g*init$Ev*init$y)))
	
	r <- init$y-X%*%beta
	V <- sum(r^2)/(init$n-ncol(A))*A	
	
	if(1-pchisq(max(beta[p]^2/V[p,p],beta[p-1]^2/V[p-1,p-1]),1) < thresh){
		mod <- lmer(update(init$formula,.~.+g+g:init$E), data = init$data,...)
		beta <- mod@fixef
		V <- vcov(mod)
		k <- pmatch(c("g","g:init$E"),names(beta))
		list(beta.main = beta[k[1]], beta.inter = beta[k[2]], 
			v.main = V[k[1],k[1]], v.inter = V[k[2],k[2]], cov.inter = V[k[1],k[2]])
		}
	
	list(beta.main = beta[p-1], beta.inter = beta[p], 
		v.main = V[p-1,p-1], v.inter = V[p,p], cov.inter = V[p,p-1])
	}

boss.fit.flm <- function(g, init, thresh=1e-7, robust=TRUE,...){
	g <- as.numeric(g)
	p <- init$p+1
	
	R <- updateR(g,init$R,init$XD)
	A <- chol2inv(R)
	
	X <- cbind(init$X,g)
	
	beta <- A%*%c(init$XY,sum(g*init$y))
	r <- init$y-X%*%beta
	if(robust){
		l <- as.vector(r)*X%*%A
		V <- t(l)%*%l	
	} else {
		V <- sum(r^2)/(init$n-ncol(A))*A	
	}
	list(beta.main=beta[p], v.main = V[p,p])
}
	
boss.fit.flm.GxE <- function(g, init, thresh=1e-7, robust=TRUE,...){
	g <- as.numeric(g)
	p <- init$p+2
	
	X <- cbind(init$X,g,g*init$E)
	R <- updateR2(cbind(g,g*init$E),init$R,init$XD)

	A <- chol2inv(R)
	
	GY <- init$y*g
	beta <- A%*%c(init$XY,sum(GY),sum(init$E*GY))
	r <- as.vector(init$y-X%*%beta)
	if(robust){
		V <- crossprod(r*X%*%A)
	} else {
		V <- sum(r^2)*A/init$n	
	}
	
	list(beta.main = beta[p-1], beta.inter = beta[p], 
		v.main = V[p-1,p-1], v.inter = V[p,p], cov.inter = V[p,p-1])
}
	
boss.fit.flr <- function(g, init, thresh = 1e-7, robust=TRUE,...){
	g <- as.numeric(g)
	p <- init$p+1
	
	#First iterate
	X <- cbind(init$X,g)
	R <- updateR(g*sqrt(init$w),R=init$R,xold=init$XD)
	A <- chol2inv(R)
	beta <- c(init$b,0)+A%*%c(init$XY,sum((init$y-init$yhat)*g))
	if(robust){
		V <- crossprod((X*as.vector(init$y-init$family$linkinv(X%*%beta)))%*%A)
	} else {
		V <- A
	}
	#The rest
	if(1-pchisq(beta[p]^2/V[p,p],1) < thresh) {
		mod <- glm.fit(X,init$y,family=binomial(), start=beta)
		beta <- coef(mod)
		if(robust){
			V <- crossprod((X*as.vector(init$y - mod$fitted))%*%chol2inv(mod$R))
		} else {
			V <- chol2inv(mod$R)	
		}
	}
	list(beta.main=beta[p], v.main=V[p,p])
}
	
boss.fit.flr.GxE <- function(g, init,thresh=1e-7, robust=TRUE,...){
	g <- as.numeric(g)
	p <- init$p+2
	
	#First iterate
	X <- cbind(init$X,g,init$E*g)
	R <- updateR2(cbind(1,init$E)*g*sqrt(init$w),R=init$R,xold=init$XD)
	A <- chol2inv(R)
	Gy <- g*(init$y-init$yhat)
	beta <- c(init$b,0,0)+A%*%c(init$XY,sum(Gy),sum(init$E*Gy))
	if(robust){
		V <- crossprod((X*as.vector(init$y-init$family$linkinv(X%*%beta)))%*%A)
	} else {
		V <- A
	}
	
	if(1-pchisq(max(beta[p]^2/V[p,p],beta[p-1]^2/V[p-1,p-1]),1) < thresh) {
		mod <- glm.fit(X,init$y,family=binomial(), start=beta)
		beta <- coef(mod)
		if(robust){
			V <- crossprod((X*as.vector(init$y - mod$fitted))%*%chol2inv(mod$R))
		} else {
			V <- chol2inv(mod$R)	
		}
	} 
	
	list(beta.main = beta[p-1], beta.inter = beta[p], 
		v.main = V[p-1,p-1], v.inter = V[p,p], cov.inter = V[p,p-1])
}

boss.fit.fgee.GxE <- function(g, init, thresh=1e-7, robust=TRUE,...){
	g <- as.numeric(g)
	p <- init$p+2
	
	Gd <- g*init$d
	Ged <- g*init$Ed
	Gy <- g*init$y
	R <- updateR2(cbind(Gd,Ged),init$R,init$XD)
	A <- chol2inv(R)
	
	X <- cbind(init$X,Gd,Ged)
	XVi <- cbind(init$XVi, g*init$v, init$Ev*g)
	beta <- A%*%c(init$XY,sum(Gy),sum(Gy*init$Ev))
	r <- as.vector(init$y-X%*%beta)
	
	if(robust){
		V <- crossprod(rowsum(XVi*r,init$id)%*%A)
	} else {
		V <- A*sum(r^2)/(init$n-init$p-1)	
	}
	
	if(1-pchisq(max(beta[p]^2/V[p,p],beta[p-1]^2/V[p-1,p-1]),1) < thresh & init$corstr != "independence"){
		m <- geese(init$y~init$X[,-1]+g+g:init$E,id=init$id, corstr=init$corstr)
		beta <- m$beta
		if(robust){
			V <- m$vbeta
		} else {
			V <- m$vbeta.naiv	
		}
	}
	list(beta.main = beta[p-1], beta.inter = beta[p], 
		v.main = V[p-1,p-1], v.inter = V[p,p], cov.inter = V[p,p-1])
}

boss.fit.fgee <- function(g, init, thresh=1e-7, robust=TRUE,...){
	g <- as.numeric(g)
	p <- init$p+1
	Gd <- g*init$dw #
	Gv <- g*init$v #
	R <- updateR(Gd,init$R,init$XD)
	A <- chol2inv(R)
	
	X <- cbind(init$X,Gd)
	XVi <- cbind(init$XVi,Gv)
	beta <- A%*%c(init$XY,sum(Gv*init$y))
	r <- as.vector(init$y-X%*%beta)
	if(robust){
		V <- crossprod((rowsum(XVi*r,init$id))%*%A)
	} else {
		V <- A*sum(r^2)/(init$n-init$p-1)	
	}
	if(1-pchisq(beta[p]^2/V[p,p],1) < thresh & init$corstr != "independence"){
		m <- geese(init$y~init$X[,-1]+g,id=init$id, corstr=init$corstr)
		beta <- m$beta
		if(robust){
			V <- m$vbeta
		} else {
			V <- m$vbeta.naiv	
		}
	}
	return(list(beta.main=beta[p], v.main = V[p,p]))
}

boss.fit.fgee.lr <- function(g, init, thresh=1e-7, robust=TRUE, ...){
	g <- as.numeric(g)
	p <- init$p+1
	
	R <- updateR(g*init$d*init$w,init$R,init$XD)
	A <- chol2inv(R)
	
	X <- cbind(init$X,g)
	XVi <- cbind(init$XVi,g*init$v*init$w)
	beta <- c(init$b,0)+A%*%c(init$XY,sum(g*init$v*init$w*(init$y-init$yhat)))
	r <- as.vector(init$y-init$yhat)	

	if(robust){	
		V <- crossprod((rowsum(XVi*r,init$id))%*%A)
	} else {
		V <- A
	}
	if(1-pchisq(beta[p]^2/V[p,p],1) < thresh){
		m <- geese(init$y~init$X[,-1]+g,id=init$id, b=beta, family=binomial(), corstr=init$corstr,...)
		beta <- m$beta
		if(robust){
			V <- m$vbeta
		} else {
			V <- m$vbeta.naiv	
		}
	}
	return(list(beta.main=beta[p], v.main = V[p,p]))
}

boss.fit.fgee.lr.GxE <- function(g, init, thresh=1e-7, robust=TRUE, ...){
	g <- as.numeric(g)
	p <- init$p+2
	
	its <- 0
	Gdw <- g*init$dw
	R <- updateR2(xnew=as.matrix(cbind(Gdw,g*init$Ed)), R=init$R, xold=as.matrix(init$XD))
	A <- chol2inv(R)
	
	X <- cbind(init$X,g,g*init$E)
	XVi <- cbind(init$XVi,g*init$vw,g*init$Ev)
	beta <- c(init$b,0,0)+A%*%t(XVi)%*%(init$y-init$yhat)
	r <- as.vector(init$y-init$yhat)	
	
	if(robust){
		V <- crossprod((rowsum(XVi*r,init$id))%*%A)
	} else {
		 V <- A
	}
	if(1-pchisq(max(beta[p]^2/V[p,p],beta[p-1]^2/V[p-1,p-1]),1) < thresh){
		m <- geese(init$y~init$X[,-1]+g+g:init$E,id=init$id, b=beta, family=binomial(), corstr=init$corstr)
		beta <- m$beta
		V <- m$vbeta
		its <- m$iterations
	}
		list(beta.main = beta[p-1], beta.inter = beta[p], 
		v.main = V[p-1,p-1], v.inter = V[p,p], cov.inter = V[p,p-1], iterations=its)
}


###wrapper for ncdf files

boss.ncdf <- function(nc, init, id.labels = NULL, g.labels = NULL, subset=NULL, gdim = 1, chunk = 1000, verbose = TRUE, outfile = NULL,...){
	if(!(gdim %in%c(1,2))) stop("Genotype dimension must be rows (1) or columns (2)")
	
	if(!is.null(g.labels)){
		g.names <- get.var.ncdf(nc = nc, varid = g.labels)
		
		if(!is.null(subset)){
			wh.g <- sort(pmatch(subset, g.names))
		} else {
			wh.g = 1:length(g.names)	
		}
	} else {
		warning("No labels given for genotype!")
	}
	if(!is.null(id.labels)){
		id <- get.var.ncdf(nc,varid = id.labels)
		g.order <- match(init$id,id)
	} else {
		warning("No ID variable given! Assuming genotype is properly ordered")
		g.order <- 1:init$n
	}
	
	nsnp <- length(g.names)
	nsnp.analyze <- length(wh.g)
	nsub <- length(id)
	nchunk <- ceiling(nsnp/chunk)
	
	if(is.null(init$E)){
		results <- matrix(NA,nrow = nsnp.analyze, ncol = 4, dimnames = list(g.names[wh.g], c("MAF","beta.main","var.main","Chi2")))
	} else {
		results <- matrix(NA,nrow = nsnp.analyze, ncol = 8, dimnames = list(g.names[wh.g], c("MAF","beta.main","beta.inter","var.main","var.inter","cov.inter","Chi2.main","Chi2.inter")))
	}
	
	i.end = 0
	for(i in 1:nchunk){
		if(i == nchunk){
			chunk <- nsnp-i.end 
		} 
		i.start = i.end + 1
		i.end = i.start + chunk-1
		
		which.snps <- (i.start:i.end) %in% wh.g
		if(gdim == 2){
			geno <- get.var.ncdf(nc,start=c(1,i.start),count=c(nsub,chunk))[,which.snps]
		} else { 
			geno <- t(get.var.ncdf(nc,start=c(i.start,1),count=c(chunk,nsub))[which.snps,])
		}
		geno.var <- apply(geno,2,var)
		geno.maf <- apply(geno,2,function(x){m <- mean(x); min(m,2-m)})
		geno <- geno[g.order,]
		if(init$method == "chol"){
			for(k in 1:chunk){
				if(is.null(init$E)){
					fit <- tryCatch(boss.fit(geno[,k],init,...), error = function(e){ cat('!'); list(beta.main = NA,v.main = NA)})
					results[i.start+k-1,] <- c(geno.maf[k], fit$beta.main, fit$v.main, fit$beta.main^2/fit$v.main)
				} else {
					fit <- tryCatch(boss.fit(geno[,k],init,...), error = function(e){ cat('!')
							 list(beta.main = NA, beta.inter = NA, v.main = NA, v.inter = NA,cov.inter = NA)})
					results[i.start+k-1,] <- c(geno.maf[k], fit$beta.main, fit$beta.inter,
						 fit$v.main, fit$v.main, fit$v.inter, fit$cov.inter, fit$beta.main^2/fit$v.main, fit$beta.inter^2/fit$v.inter)
				}
			}
		} else {
			fit <- boss.fit(geno,init)
			results[i.start:i.end,] <- cbind(geno.maf,fit$beta.main,fit$v.main, fit$chi2)
		}
		if(verbose){
			cat(".")
			if(i%%10 == 0) cat(i.end,"\n")
		}	
	}
	if(is.null(outfile)){
		return(results)	
	} else {
		write.csv(results,file = paste(outfile,".csv"))
	}
}

