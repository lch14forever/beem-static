######### internal functions #########

##' @title tss
##'
##' @param x count matrix (each OTU in one row)
##' @description Calculate the relative abundances by dividing the total abundance (Total Sum Sacling)
##' @author Chenhao Li, Niranjan Nagarajan
tss <- function(x){ apply(x, 2, function(x) x/sum(x)) }


##' @title css
##'
##' @param m matrix of data (variables in columns, measurements in rows)
##' @param p quantile used for normalization (default: 0.5)
##' @description Function to perform cumulative sum scaling (CSS)
##' @author Chenhao Li, Niranjan Nagarajan
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
##' @param alpha the alpha parameter for elastic net (1:lasso [default], 0:ridge)
##' @param lambda.choice 1: use lambda.1se for analysis, 2: use lambda.min for analysis, 3: mixed
##' @import glmnet
##' @description Infer parameters with blasso
##' @author Chenhao Li, Niranjan Nagarajan
infer <- function(Y, X, method='glmnet', intercept=FALSE, seed=0, alpha=1, lambda.choice=1){
    set.seed(seed)
    if(method=='lm'){
        return(as.numeric(coef(lm(Y~X+0))))
    }
    if(method=='glmnet'){
        lambda.init <- rev(10^(seq(-9,-1, length.out = 20)))
        penalty <- c(0,rep(1, ncol(X)-1))
        maxmin <- median(Y) + 5 * IQR(Y) %*% c(1,-1)
        idx <- Y <= maxmin[1] & Y >=maxmin[2]
        fit <- cv.glmnet(X[idx,], Y[idx], intercept=intercept, lambda=lambda.init,
                         penalty.factor=penalty, alpha=alpha)

        lambda <- rev(exp(seq(log(fit$lambda.1se/20), log(fit$lambda.1se*20), length.out = 100)))

        fit <- cv.glmnet(X[idx,], Y[idx], intercept=intercept, lambda=lambda,
                         penalty.factor=penalty, alpha=alpha)
        if(lambda.choice == 1){
            s = 'lambda.1se'
        }else if(lambda.choice == 2){
            s = 'lambda.min'
        }else{
            s = (1-lambda.choice)*fit$lambda.min + lambda.choice*fit$lambda.1se
        }
        coefs <- coef(fit, s=s)[-1]

        e2 <- as.numeric((Y-(X %*% coefs)[,1])^2)
        return(c(coefs, e2))
    }
}

##' @title norm
##'
##' @param Y response
##' @param X predictors
##' @description estimate biomass with linear regression
##' @author Chenhao Li, Niranjan Nagarajan
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
##' @author Chenhao Li, Niranjan Nagarajan
norm <- function(a, b, x){
    res <- -a / (b %*% x)
    res <- res[x!=0]
    if(all(res < 0)){
        m <- abs(max(res))
    } else{
        m <- median(res[res>0])
    }
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
##' @importFrom doMC registerDoMC
##' @import foreach
##' @description E-part of BEEM, estimate model parameters with inferred m
##' @author Chenhao Li, Niranjan Nagarajan
func.E <- function(dat.tss, m, sample.filter, ncpu=4, center=FALSE, ...){
    ## infer parameter for each OTU
    registerDoMC(ncpu)
    p <- nrow(dat.tss)

    res <- foreach(i=1:p, .combine=rbind) %dopar% {
        fil <- dat.tss[i,]!=0 & !sample.filter[i,]
        X <- t(rbind(1/m, dat.tss[-i,])[,fil])
        Y <- dat.tss[i, fil]
        if(center){
            Y <- Y-mean(Y)
            X <- X-rowMeans(X)
        }
        theta <- rep(0, p+1)
        theta[i+1] <- -1 ## -beta_{ii}/beta_{ii}
        tmp <- infer(Y,X,...)
        theta[-(i+1)] <- tmp[1:p]
        e2 <- rep(NA, ncol(dat.tss))
        e2[fil] <- tmp[-(1:p)]
        c(theta, e2)
    }
    list(a=res[,1], b=postProcess(res[,1:p+1]), e2=res[,-(1:(p+1))])
}


##' @title func.M
##'
##' @param dat.tss relative abundance matrix (each OTU in one row)
##' @param a estimated growth rate values (scaled by self-interaction)
##' @param b estimated interaction matrix (scaled by self-interaction)
##' @param ncpu number of CPUs (default: 4)
##' @importFrom doMC registerDoMC
##' @import foreach
##' @description M-part of BEEM, estimate biomass with inferred model parameters
##' @author Chenhao Li, Niranjan Nagarajan
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
##' @author Chenhao Li, Niranjan Nagarajan
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
##' @author Chenhao Li, Niranjan Nagarajan
postProcess <- function(b, dev=1e-5){
    b[abs(b)<dev] <- 0
    b
}


##' @title formatOutput
##' @param a scaled growth rates
##' @param b scaled interaction matrix
##' @param vnames variable names
##' @description Function to convert parameter vector a and matrix b to MDSINE's output format
##' @author Chenhao Li, Niranjan Nagarajan
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

##' @title detectBadSamples
##'
##' @param err the errors estimated from regression
##' @param threshold threshold to filter out samples
##' @author Chenhao Li, Niranjan Nagarajan
detectBadSamples <- function(err, threshold){
    score <- abs(err - median(err, na.rm = TRUE))/IQR(err, na.rm = TRUE)
    score[is.na(score)] <- Inf
    return(score > threshold)
}

######### Exported functions #########

##' @title beem2param
##'
##' @param beem a BEEM object
##' @description extract parameter estimates from a BEEM object
##' @export
##' @author Chenhao Li, Niranjan Nagarajan
beem2param <- function(beem){
    p <- ncol(beem$err.m)
    tmp <- beem$trace.p[, ncol(beem$trace.p)]
    a.est <- tmp[1:p]
    b.est <- matrix(tmp[-c(1:p)], p, p)
    return (list(a.est=a.est, b.est=b.est))
}

##' @title beem2biomass
##'
##' @param beem a BEEM object
##' @description extract biomass estimates from a BEEM object
##' @export
##' @author Chenhao Li, Niranjan Nagarajan
beem2biomass <- function(beem){
    niter <- ncol(beem$trace.m)
    return (beem$trace.m[,niter])
}


##' @title func.EM
##'
##' @param dat OTU count/relative abundance matrix (each OTU in one row)
##' @param ncpu number of CPUs (default: 4)
##' @param scaling a scaling factor to keep the median of all biomass constant (default: 1000)
##' @param dev deviation of the error (for one sample) from the model to be excluded (default: Inf - all the samples will be considered)
##' @param max.iter maximal number of iterations (default 30)
##' @param warm.iter number of iterations to run before removing any samples (default: run until convergence and start to remove samples)
##' @param lambda.choice 1: use lambda.1se for LASSO, 2: use lambda.min for LASSO, a number between (0, 1): this will select a lambda according to (1-lambda.choice)*lambda.min + lambda.choice*lambda.1se
##' @param alpha The alpha parameter for the Elastic Net model (1-LASSO [default], 0-RIDGE)
##' @param refresh.iter refresh the removed samples every X iterations (default: 3)
##' @param debug output debugging information (default FALSE)
##' @description Iteratively estimating scaled parameters and biomass
##' @export
##' @author Chenhao Li, Niranjan Nagarajan
func.EM <- function(dat, ncpu=4, scaling=10000, dev=Inf, max.iter=30,
                    warm.iter=NULL, lambda.choice=1, alpha=1, debug=FALSE, refresh.iter=3){

    ## pre-processing
    dat.init <- preProcess(dat, dev=0)
    dat.tss <- dat.init$tss
    spNames <- rownames(dat)

    ## ensure valid samples
    temp <- colSums(dat.tss, na.rm = TRUE) == 0
    if(any(temp)){
        stop(paste0('Sample ', which(temp), ' has zero total abudance...'))
    }

    ## initialization
    sample.filter.iter <- dat.init$sample.filter
    tmp <- css(t(dat.tss))$normFactors
    m.iter <- scaling * tmp/median(tmp)
    ## m.iter <- rnorm(length(m.iter), scaling, scaling/10) ## start with random biomass
    trace.m <- matrix(m.iter)
    trace.p <- apply(expand.grid(spNames, spNames), 1, function(x) paste0(x[2], '->', x[1]))
    trace.p <- c(paste0(NA, '->',spNames), trace.p)
    trace.p <- data.frame(name=trace.p)

    ## flags
    remove_non_eq <- FALSE
    removeIter <- 0

    ## constants
    m1 <- matrix(rep(1,nrow(sample.filter.iter)))
    ## EM
    for(iter in 1:max.iter){
        message(paste0("############# Run for iteration ", iter,": #############"))
        message("E-step: estimating scaled parameters...")
        ##if(iter==1) method <- 'lm'
        if(iter>=1) method <- 'glmnet'
        ## TODO: is a switch from lasso to elastic net needed?
        tmp.p <- func.E(dat.tss, m.iter, sample.filter.iter, ncpu, method=method, lambda.choice=lambda.choice, alpha=alpha)
        err.p <- tmp.p$e2
        if(debug){
            print(rowSums(tmp.p$b!=0))
            plot(1-rowSums(err.p, na.rm = TRUE)/apply(dat.tss, 1, function(x) sum((x[x!=0]-mean(x[x!=0]))^2)))
            abline(h=0.5)
        }
        message("M-step: estimating biomass...")
        tmp.m <- func.M(dat.tss, tmp.p$a, tmp.p$b, ncpu)
        m.iter <- tmp.m[,1]
        err.m <- tmp.m[,-1]

        if(remove_non_eq){
            ## clear up removed samples every X iterations
            if (iter %% refresh.iter == 0 ) {
                sample.filter.iter <- dat.init$sample.filter
            }
            bad.samples <- detectBadSamples(apply(err.p, 2, median, na.rm = TRUE), dev)
            sample.filter.iter <- (m1 %*% bad.samples)>0  | sample.filter.iter
            message(paste0("Number of samples removed (detected to be non-static): ",
                           sum(colSums(sample.filter.iter)>0)))
            removeIter <- removeIter + 1
        }

        m.iter <- m.iter*scaling/median(m.iter[colSums(sample.filter.iter)==0], na.rm=TRUE)
        trace.m <- cbind(trace.m, m.iter)
        trace.p <- cbind(trace.p, formatOutput(tmp.p$a, tmp.p$b, spNames)$value)
        criterion <- median(abs((trace.m[, iter+1] - trace.m[, iter])/trace.m[, iter]))<1e-3
        if (!is.null(warm.iter) && iter > warm.iter && !remove_non_eq){
            message("Start to detect and remove bad samples...")
            remove_non_eq <- TRUE
        }
        if (iter > 5 && !remove_non_eq && criterion && is.finite(dev)) {
            message("Converged and start to detect and remove bad samples...")
            remove_non_eq <- TRUE
        }
        if (((removeIter > 5 && remove_non_eq) || is.infinite(dev)) && criterion) {
            message("Converged!")
            break
        }
    }
    list(trace.m=trace.m, trace.p=trace.p, err.m=err.m, err.p=err.p,
         sample2rm = which(colSums(sample.filter.iter) > 0))
}

