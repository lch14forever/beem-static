require('foreach')
require('doMC')
suppressMessages(library(glmnet))

######### internal functions #########


##' @title tss
##' 
##' @param x count matrix (each OTU in one row)
##' @description Calculate the relative abundances by dividing the total abundance (Total Sum Sacling)
tss <- function(x){ apply(x, 2, function(x) x/sum(x)) }


##' @title css
##' 
##' @param m matrix of data (variables in columns, measurements in rows)
##' @param p quantile used for normalization (default: 0.5)
##' @description Function to perform cumulative sum scaling (CSS)
css <- function(m, p=0.5){
    m[m==0] <- NA
    ## find the quantile in each sample
    quant <- apply(m,1,function(x) quantile(x, p=p, na.rm = TRUE))
    ## calculate normalization factor
    f <- rowSums(m*sweep(m, 1, quant, '<='),na.rm = TRUE)
    nf <- f/exp(mean(log(f)))
    dat.css <- sweep(m,1,nf, '/')
    return(list("normCounts" = dat.css, "normFactors"=nf))
}


##' @title infer
##'
##' @param Y response
##' @param X predictors
##' @param method blasso or lm
##' @param intercept whether to include the intercept
##' @param seed seed
##' @description Infer parameters with blasso
infer <- function(Y, X, method='glmnet', intercept=FALSE, seed=0){
    ##set.seed(seed)
    ## if(method=='monomvn'){
    ##     res <- tryCatch({
    ##         bl.fit <- blasso(X, Y, verb = 0, T=1000)
    ##         b <- apply(bl.fit$beta[-c(1:500),],2, median)
    ##         return(1)
    ##     },
    ##     error=function(x){
    ##         message("blasso error caught, using glmnet lasso instead")
    ##         return(NA)
    ##     })
    ##     if(is.na(res)){method <- 'glmnet'}else{return(b)}            
    ## }
    ## if(method=='ebglmnet'){
    ##     library(EBglmnet)
    ##     tmp1 <- capture.output(tmp <- cv.EBglmnet(X,Y, prior='lasso'))
    ##     theta <- rep(0, ncol(X))
    ##     theta[tmp$fit[,1]] <- tmp$fit[,3]
    ##     return(theta)
    ## }
    if(method=='lm'){
        return(as.numeric(coef(lm(Y~X+0))))
    }
    if(method=='glmnet'){
        lambda.init <- 10^(seq(-7,-1))
        penalty <- c(0,rep(1, ncol(X)-1))
        ## maxmin <- median(Y) + 2* IQR(Y) %*% c(1,-1)
        ## idx <- Y <= maxmin[1] & Y >=maxmin[2]
        fit <- cv.glmnet(X, Y, intercept=intercept, lambda=lambda.init,
                         penalty.factor=penalty)
        lambda <- seq(fit$lambda.1se/10, fit$lambda.1se*10, fit$lambda.1se/10)
        fit <- cv.glmnet(X, Y, intercept=intercept, lambda=lambda,
                         penalty.factor=penalty)        
        return(as.numeric(coef(fit))[-1])
    }
}

##' @title norm
##'
##' @param Y response
##' @param X predictors
##' @description estimate biomass with linear regression
## norm <- function(Y, X){
##     wrong.sign <- sign(X) != sign(Y)
##     outliers <- abs(Y - median(Y)) > 2 * IQR(Y) | abs(X - median(X)) > 2 * IQR(X) | wrong.sign
##     if(all(outliers)) outliers <- wrong.sign
##     if(all(outliers)) return(median(abs(Y/X)))
##     model <- lm(Y[!outliers]~X[!outliers]-1)
##     abs(model$coeff[1])
##     ##median(Y/X)
## }

##' @title norm
##'
##' @param a estimated growth rate values (scaled by self-interaction)
##' @param b estimated interaction matrix (scaled by self-interaction)
##' @param x abundances in one sample
##' @description estimate biomass with linear regression
norm <- function(a, b, x){
    res <- -a / (b %*% x)
    m <- median(abs(res[x!=0]))
    err <- (m * (b %*% x) + a)[,1]/a
    err[x==0] <- 0
    c(m, err)
}


##' @title func.E
##'
##' @param dat.tss relative abundances matrix (each OTU in one row)
##' @param sample.filter filter out samples contain outliers in Y 
##' @param m estimated biomass values
##' @param ncpu number of CPUs (default: 4)
##' @param center center data or not
##' @description E-part of BEEM, estimate model parameters with inferred m
func.E <- function(dat.tss, m, sample.filter, ncpu=4, center=FALSE, ...){
    ## infer parameter for each OTU
    registerDoMC(ncpu)
    res <- foreach(i=1:nrow(dat.tss), .combine=rbind) %dopar% {
        fil <- dat.tss[i,]!=0 & !sample.filter[i,]
        X <- t(rbind(1/m, dat.tss[-i,])[,fil])
        Y <- dat.tss[i, fil]
        if(center){
            Y <- Y-mean(Y)
            X <- X-rowMeans(X)
        }
        theta <- rep(0, nrow(dat.tss)+1)
        theta[i+1] <- -1 ## -beta_{ii}/beta_{ii}
        theta[-(i+1)] <- infer(Y, X, ...)
        theta
    }    
    list(a=res[,1], b=postProcess(res[,-1]))
}


##' @title func.M
##'
##' @param dat.tss relative abundance matrix (each OTU in one row)
##' @param a estimated growth rate values (scaled by self-interaction)
##' @param b estimated interaction matrix (scaled by self-interaction)
##' @param ncpu number of CPUs (default: 4)
##' @description M-part of BEEM, estimate biomass with inferred model parameters
func.M <- function(dat.tss, a, b, ncpu=4,...){
    registerDoMC(ncpu)    
    foreach(i=1:ncol(dat.tss), .combine=rbind) %dopar%{
        x <- dat.tss[,i]
        norm(a, b, x)
    }
}

##' @title preProcess
##'
##' @param dat OTU count/relative abundance matrix (each OTU in one row)
##' @param dev dev * IQR from median will be filtered out (default: Inf, nothing to remove)
##' @param ncpu number of CPUs (default: 4)
##' @description pre-process data
preProcess <- function(dat, dev=1){
    ## filter out species abundances that are too low
    detection_limit <- 1e-4
    dat <- tss(dat)
    dat[dat<detection_limit] <- 0    
    ##rowMed <- apply(dat, 1, function(x) median(x[x!=0]))
    ##rowIQR <- apply(dat, 1, function(x) IQR(x[x!=0]))
    rowMAD <- apply(dat, 1, function(x) mad(x[x!=0]))
    ## too close to zero
    sample.filter <- dat/rowMAD < dev
    return(list(tss=dat, sample.filter=sample.filter))
}

##' @title postProcess
##'
##' @param b scaled interaction matrix (diagonal/self-interaction is -1)
##' @param dev scale of the interaction < dev * 1, the value will be set to 0
##' @description post-process parameters to remove small entries
postProcess <- function(b, dev=1e-5){
    b[abs(b)<dev] <- 0
    b
}



##' @param a
##' @param b
##' @param vnames variable names
##' @description Function to convert parameter vector a and matrix b to MDSINE's output format
formatOutput <- function(a, b, vnames){
    suppressMessages(require(reshape2))
    p <- length(a)
    param <- data.frame(parameter_type=c(rep("growth_rate",p), rep("interaction",p*p)))
    param$source_taxon <- c(rep(NA,p),rep(vnames, each = p))
    param$target_taxon <- vnames
    param$value <- 0
    param$significance <- NA
    tmp.b <- melt(t(b))
    tmp.b <- tmp.b[order(tmp.b$Var1),]
    param$value <- c(a, tmp.b$value)
    param
}

##' @title func.EM
##' 
##' @param dat OTU count/relative abundance matrix (each OTU in one row)
##' @param ncpu number of CPUs (default: 4)
##' @description Iteratively estimating scaled parameters and biomass
##' @export
func.EM <- function(dat, ncpu=4, scaling=10000, dev=2, max.iter=30, refine.start.iter=max.iter/2){
    ## pre-processing    
    tmp <- preProcess(dat, dev=0)
    dat.tss <- tmp$tss
    spNames <- rownames(dat)
    
    ## initialization
    sample.filter.iter <- tmp$sample.filter
    tmp <- css(t(dat.tss))$normFactors
    m.iter <- scaling * tmp/median(tmp) 
    trace.m <- matrix(m.iter)
    trace.p <- matrix(,nrow=nrow(dat)^2+nrow(dat))
    trace.p <- apply(expand.grid(spNames, spNames), 1, function(x) paste0(x[2], '->', x[1]))
    trace.p <- c(paste0(NA, '->',spNames), trace.p)
    trace.p <- data.frame(name=trace.p)
    trace.e.pc1 <- matrix(nrow=length(m.iter), ncol=0)
    trace.e.pc2 <- matrix(nrow=length(m.iter), ncol=0)
    ## EM
    for(iter in 1:max.iter){
        message(paste0("############# Run for iteration ", iter,": #############"))
        message("E-step: estimating scaled parameters...")
        ##if(iter==1) method <- 'lm'
        if(iter>=1) method <- 'glmnet'
        tmp.p <- func.E(dat.tss, m.iter, sample.filter.iter, ncpu, method=method)
        message("M-step: estimating biomass...")
        tmp.m <- func.M(dat.tss, tmp.p$a, tmp.p$b, ncpu)
        m.iter <- tmp.m[,1]
        err <- tmp.m[,-1]
        pca <- prcomp( log(abs(err) + 1) )
        trace.e.pc1 <- cbind(trace.e.pc1, pca$x[,1])
        trace.e.pc2 <- cbind(trace.e.pc2, pca$x[,1])
        ##err <- err/sd(err[err!=0])
        if(iter > refine.start.iter){
            bad.samples <- rowSums(apply(pca$x[,1:2], 2, function(x) abs(x-median(x))/mad(x) > dev)) > 0
            sample.filter.iter <- (matrix(rep(1,nrow(sample.filter.iter))) %*% bad.samples)>0  | sample.filter.iter
            ##sample.filter.iter <- t(abs(err)) > dev | sample.filter.iter
            
            message(paste0("Number of samples removed (detected to be non-static): ",
                                   sum(colSums(sample.filter.iter)>0)))
        }
        m.iter <- m.iter*scaling/median(m.iter[colSums(sample.filter.iter)==0])
        trace.m <- cbind(trace.m, m.iter)
        trace.p <- cbind(trace.p, formatOutput(tmp.p$a, tmp.p$b, spNames)$value)
    }
    list(trace.m=trace.m, trace.p=trace.p, err=err,
         sample2rm = which(colSums(sample.filter.iter) > 0))
}

