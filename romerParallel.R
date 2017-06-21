# Gabriel Hoffman
# April 3, 2017
#
# Parallelize romer when statistic="mean"

library(compiler)
library(foreach)
library(limma)

##	ROMER.R
romerParallel <- function(y,...) UseMethod("romerParallel")

romerParallel.default <- cmpfun(function(y,index,design=NULL,contrast=ncol(design),array.weights=NULL,block=NULL,correlation=NULL,set.statistic="mean",nrot=9999,shrink.resid=TRUE,n_slices = 100,...)
#	rotation mean-rank version of GSEA (gene set enrichment analysis) for linear models
#	Gordon Smyth and Yifang Hu
#	27 March 2009.	Last modified 3 May 2015.
{
#	Issue warning if extra arguments found
	dots <- names(list(...))
	if(length(dots)) warning("Extra arguments disregarded: ",sQuote(dots))

#	Check y
	y <- as.matrix(y)
	ngenes <- nrow(y)
	n <- ncol(y)

#	Check index
	if(!is.list(index)) index <- list(set=index)
	nsets <- length(index)
	if(nsets==0) stop("index is empty")
	SetSizes <- unlist(lapply(index,length))

#	Check design
	if(is.null(design))
		stop("design matrix not specified")
	else {
		design <- as.matrix(design)
		if(mode(design) != "numeric") stop("design must be a numeric matrix")
	}
	if(nrow(design) != n) stop("row dimension of design matrix must match column dimension of data")
	ne <- nonEstimable(design)
	if(!is.null(ne)) cat("Coefficients not estimable:",paste(ne,collapse=" "),"\n")
	p <- ncol(design)
	if(p<2L) stop("design needs at least two columns")
	p0 <- p-1L
	d <- n-p

#	Reform design matrix so that contrast of interest is last column
	if(length(contrast) == 1L) {
		contrast <- round(contrast)
		if(contrast < p) {
			i <- 1L:p
			design <- design[,c(i[-contrast],contrast)]
		}
	} else {
		design <- contrastAsCoef(design=design, contrast=contrast, first=FALSE)$design
	}

#	Divide out array weights, if they exist
	if(!is.null(array.weights)) {
		if(any(array.weights <= 0)) stop("array.weights must be positive")
		if(length(array.weights) != n) stop("Length of array.weights doesn't match number of array")
		design <- design*sqrt(array.weights)
		y <- t(t(y)*sqrt(array.weights))
	}

#	Divide out correlation structure, it is exists
	if(!is.null(block)) {
		block <- as.vector(block)
		if (length(block) != n) stop("Length of block does not match number of arrays")
		ub <- unique(block)
		nblocks <- length(ub)
		Z <- matrix(block,n,nblocks) == matrix(ub,n,nblocks,byrow = TRUE)
		cormatrix <- Z %*% (correlation * t(Z))
		diag(cormatrix) <- 1
		R <- chol(cormatrix)
		y <- t(backsolve(R, t(y), transpose = TRUE))
		design <- backsolve(R, design, transpose = TRUE)
	}

# 	Fit linear model to all genes
	qr <- qr(design)
	signc <- sign(qr$qr[p,p])
	effects <- qr.qty(qr,t(y))

#	Sample variances
	s2 <- colMeans(effects[-(1:p),,drop=FALSE]^2)

# 	Estimate global hyper-parameters s0 and d0
	sv <- squeezeVar(s2,df=d)
	d0 <- sv$df.prior
	s02 <- sv$var.prior
	sd.post <- sqrt(sv$var.post)

#	t-statistics and effects (orthogonal residuals)
	Y <- effects[-(1:p0),,drop=FALSE]
	YY <- colSums(Y^2)
	B <- Y[1,]
	modt <- signc*B/sd.post

	##	EMP BAYES SHRINK FIRST EFFECT TO REMOVE SYSTEMATIC COMPONENT
	if(shrink.resid) {

	#	Estimate hyperparameter p0 (proportion of DE genes)
		p.value <- 2*pt(abs(modt),df=d0+d,lower.tail=FALSE)
		proportion <- 1-propTrueNull(p.value) # proportion of DE probes

	#	Estimate hyperparameter v0 (var.prior, variance of true logFC)
	 	stdev.unscaled <- rep_len(1/abs(qr$qr[qr$rank,qr$rank]),ngenes)
	 	var.unscaled <- stdev.unscaled^2
		df.total <- rep_len(d,ngenes) + sv$df.prior
		stdev.coef.lim <- c(0.1, 4)
		var.prior.lim <- stdev.coef.lim^2/sv$var.prior
		var.prior <- tmixture.vector(modt, stdev.unscaled, df.total, proportion, var.prior.lim)
		if (any(is.na(var.prior))) {
			var.prior[is.na(var.prior)] <- 1/sv$var.prior
			warning("Estimation of var.prior failed - set to default value")
		}

	#	Estimate posterior probability of DE
		r <- (var.unscaled + var.prior)/var.unscaled
		if (sv$df.prior > 10^6)
			kernel <- modt^2 * (1 - 1/r)/2
		else
			kernel <- (1 + df.total)/2 * log((modt^2 + df.total)/(modt^2/r + df.total))
		lods <- log(proportion/(1 - proportion)) - log(r)/2 + kernel
		ProbDE <- exp(lods)/(1+exp(lods))

	#	Shrink contrast to be like a residual
		Y[1,] <- Y[1,]*sqrt(var.unscaled/(var.unscaled+var.prior*ProbDE))

	}
##	END SHRINK

	set.statistic <- match.arg(set.statistic,c("mean","floormean","mean50"))

	if(set.statistic=="mean") {

		#	Observed rankings for each set
		obs.ranks <- matrix(0,ngenes,3)
		obs.ranks[,1] <- rank(modt)
		obs.ranks[,2] <- ngenes-obs.ranks[,1]+1
		obs.ranks[,3] <- rank(abs(modt))

		AllIndices <- unlist(index)
		Set <- rep(1:nsets,SetSizes)
		obs.set.ranks <- rowsum(obs.ranks[AllIndices,],group=Set,reorder=FALSE)
		obs.set.ranks <- obs.set.ranks / SetSizes

		#	Random rotations to simulate null hypothesis
		rot.ranks <- obs.ranks
		# p.value <- matrix(0,nrow=nsets,ncol=3)

		# n_slices = 100

		p.value = foreach( j = 1:n_slices, .combine=function(A, B){ A+B}) %dopar% {
			cat("Slice:", j, "\n")
			foreach(i=1:ceiling(nrot/n_slices), .combine=function(A, B){ A+B}) %do% {
				R <- matrix(rnorm((d+1)),1,d+1)
				R <- R/sqrt(rowSums(R^2))
				Br <- R %*% Y
				s2r <- (YY-Br^2)/d

				if(is.finite(d0)){
					sdr.post <- sqrt((d0*s02+d*s2r)/(d0+d))
				}else{
					sdr.post <- sqrt(s02)
				}

				modtr <- signc*Br/sdr.post
			
				rot.ranks[,1] <- rank(modtr)
				rot.ranks[,2] <- ngenes-rot.ranks[,1]+1
				rot.ranks[,3] <- rank(abs(modtr))

				rot.set.ranks <- rowsum(rot.ranks[AllIndices,],group=Set,reorder=FALSE)
				rot.set.ranks <- rot.set.ranks / SetSizes
				# p.value <- p.value + (rot.set.ranks >= obs.set.ranks)

				rm(R, Br, s2r, sdr.post, modtr)#, rot.ranks)
				(rot.set.ranks >= obs.set.ranks)
			}
		}	

	} # end "mean"

	if(set.statistic=="floormean") {

	#	Observed rankings for each set
		obs.ranks <- matrix(0,ngenes,3)
		obs.ranks[,1] <- rank(pmax(modt,0))
		obs.ranks[,2] <- rank(pmax(-modt,0))
		obs.ranks[,3] <- rank(pmax(abs(modt),1))

		AllIndices <- unlist(index)
		Set <- rep(1:nsets,SetSizes)
		obs.set.ranks <- rowsum(obs.ranks[AllIndices,],group=Set,reorder=FALSE)
		obs.set.ranks <- obs.set.ranks / SetSizes

	#	Random rotations to simulate null hypothesis
		rot.ranks <- obs.ranks
		p.value <- matrix(0,nrow=nsets,ncol=3)
		for(i in 1:nrot)
		{
			R <- matrix(rnorm((d+1)),1,d+1)
			R <- R/sqrt(rowSums(R^2))
			Br <- R %*% Y
			s2r <- (YY-Br^2)/d

			if(is.finite(d0))
				sdr.post <- sqrt((d0*s02+d*s2r)/(d0+d))
			else
				sdr.post <- sqrt(s02)

			modtr <- signc*Br/sdr.post
		
			rot.ranks[,1] <- rank(pmax(modtr,0))
			rot.ranks[,2] <- rank(pmax(-modtr,0))
			rot.ranks[,3] <- rank(pmax(abs(modtr),1))

			rot.set.ranks <- rowsum(rot.ranks[AllIndices,],group=Set,reorder=FALSE)
			rot.set.ranks <- rot.set.ranks / SetSizes
			p.value <- p.value + (rot.set.ranks >= obs.set.ranks)
		}	
	} # end "floormean"

	if(set.statistic=="mean50") {
	
	#	Observed rankings for each set
		s.r <- rank(modt)
		s.abs.r <- rank(abs(modt))

		s.rank.mixed <- rep(0,nsets)
		s.rank.up <- rep(0,nsets)
		s.rank.down <- rep(0,nsets)

		m <- floor((SetSizes+1)/2)
		for(i in 1:nsets)
		{
			mh <- limma:::.meanHalf(s.r[index[[i]]],m[i])
			s.rank.up[i] <- mh[2]	
			s.rank.down[i] <- mh[1]
			s.rank.mixed[i] <- limma:::.meanHalf(s.abs.r[index[[i]]],m[i])[2]
		}	

	#	Random rotations
		p.value <- matrix(rep(0,nsets*3),nrow=nsets,ncol=3)
		for(i in 1:nrot)
		{
			R <- matrix(rnorm((d+1)),1,d+1)
			R <- R/sqrt(rowSums(R^2))
			Br <- R %*% Y
			s2r <- (YY-Br^2)/d

			if(is.finite(d0))
				sdr.post <- sqrt((d0*s02+d*s2r)/(d0+d))
			else
				sdr.post <- sqrt(s02)

			modtr <- signc*Br/sdr.post

			s.r.2 <- rank(modtr)
			s.abs.r.2 <- rank(abs(modtr))
		
			for(j in 1:nsets)
			{
				mh.2 <- limma:::.meanHalf(s.r.2[index[[j]]],m[j])

				s.rank.up.2 <- mh.2[2]
				s.rank.down.2 <- mh.2[1]
				s.rank.mixed.2 <- limma:::.meanHalf(s.abs.r.2[index[[j]]],m[j])[2]
			
				if(s.rank.up.2 >= s.rank.up[j]) p.value[j,1] <- p.value[j,1]+1
				if(s.rank.down.2 <= s.rank.down[j]) p.value[j,2] <- p.value[j,2]+1
				if(s.rank.mixed.2 >= s.rank.mixed[j]) p.value[j,3] <- p.value[j,3]+1
			}
		}	
	} # end "mean50"

	p.value <- (p.value+1)/(nrot+1)
	colnames(p.value) <- c("Up","Down","Mixed")
	SetNames <- names(index)
	if(is.null(SetNames))
		rownames(p.value) <- 1:nsets
	else
		rownames(p.value) <- SetNames
	cbind(NGenes=SetSizes,p.value)
})

