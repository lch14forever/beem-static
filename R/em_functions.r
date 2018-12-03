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
##' @param lambda.init user provides initial lambda values
##' @param lambda.choice 1: use lambda.1se for LASSO, 2: use lambda.min for LASSO, a number between (0, 1): this will select a lambda according to (1-lambda.choice)*lambda.min + lambda.choice*lambda.1se
##' @param nfolds number of folds for glmnet cross-validation
##' @import glmnet
##' @description Infer parameters with blasso
##' @author Chenhao Li, Niranjan Nagarajan
infer <- function(Y, X, method='glmnet', intercept=FALSE, seed=0, alpha=1, lambda.init=NULL, lambda.choice=1, nfolds=5){
    set.seed(seed)
    if(method=='lm'){
        return(as.numeric(coef(lm(Y~X+0))))
    }
    if(method=='glmnet'){
        penalty <- c(0,rep(1, ncol(X)-1))
        maxmin <- median(Y) + 5 * IQR(Y) %*% c(1,-1)
        idx <- Y <= maxmin[1] & Y >=maxmin[2]
        if(is.null(lambda.init)){
            ## select a lambda range first
            lambda.init <- rev(10^(seq(-9,-1, length.out = 20)))
            fit <- cv.glmnet(X[idx,], Y[idx], intercept=intercept, lambda=lambda.init, nfolds=nfolds,
                         penalty.factor=penalty, alpha=alpha)
            lambda <- rev(exp(seq(log(fit$lambda.min/10), log(fit$lambda.min*10), length.out = 50)))
        }else{
            lambda <- rev(exp(seq(log(lambda.init/2), log(lambda.init), length.out = 10)))
        }

        # ## Iteratively fit models
        # lambdas <- NULL #initialize
        # for (i in 1:5) {
        #     fit <- cv.glmnet(X[idx,], Y[idx], intercept=intercept, lambda=lambda,  nfolds=nfolds,
        #                      penalty.factor=penalty, alpha=alpha)
        #     if(lambda.choice == 1){
        #         s = fit$lambda.1se
        #     }else if(lambda.choice == 2){
        #         s = fit$lambda.min
        #     }else{
        #         s = (1-lambda.choice)*fit$lambda.min + lambda.choice*fit$lambda.1se
        #     }
        #
        #     df <- data.frame(s, median(fit$cvm) ) #can use median for CVM also
        #     lambdas <- rbind(lambdas, df)
        # }
        # ## Select best lambda:
        # s <- median(lambdas[, 1])

        fit <- cv.glmnet(X[idx,], Y[idx], intercept=intercept, lambda=lambda,  nfolds=nfolds,
                         penalty.factor=penalty, alpha=alpha)
        if(lambda.choice == 1){
            s = fit$lambda.1se
        }else if(lambda.choice == 2){
            s = fit$lambda.min
        }else{
            s = (1-lambda.choice)*fit$lambda.min + lambda.choice*fit$lambda.1se
        }
        coefs <- coef(fit, s=s)[-1]
        e2 <- as.numeric((Y-(X %*% coefs)[,1])^2)
        rm(.Random.seed, envir=.GlobalEnv)
        return(c(coefs, e2, s))
    }
}

# ##' @title boot_stat
# ##'
# ##' @param data data in the format of cbind(Y, X)
# ##' @param indices indices for bootstrapping
# ##' @description bootstrapping the inference process
# ##' @author Chenhao Li, Niranjan Nagarajan
# bs <- function(data, indices) {
#     X_s <- data[indices, -1] # allows boot to select sample
#     Y_s <- data[indices, 1]
#     fit <- infer(Y_s, X_s)
#     return(fit[1:ncol(X_s)])
# }
# ##' @title func.E.boot
# ##'
# ##' @param dat.tss relative abundances matrix (each OTU in one row)
# ##' @param sample.filter filter out samples contain outliers in Y
# ##' @param m estimated biomass values
# ##' @param ncpu number of CPUs (default: 4)
# ##' @importFrom boot boot
# ##' @importFrom doMC registerDoMC
# ##' @import foreach
# ##' @description bootstrapped E-step
# ##' @author Chenhao Li, Niranjan Nagarajan
# func.E.boot <- function(dat.tss, m, sample.filter, ncpu=4, ...){
#     registerDoMC(ncpu)
#     p <- nrow(dat.tss)
#     res <- foreach(i=1:p, .combine=rbind) %do% {
#         message(paste0("Bootstrapping for species ", i))
#         fil <- dat.tss[i,]!=0 & !sample.filter[i,]
#         X <- t(rbind(1/m, dat.tss[-i,])[,fil])
#         Y <- dat.tss[i, fil]
#         theta <- rep(0, p+1)
#         stab <- rep(1, p)
#         theta[i+1] <- -1 ## -beta_{ii}/beta_{ii}
#         tmp <- boot(data=cbind(Y,X), statistic=bs, R=100)
#         theta[-(i+1)] <- apply(tmp$t, 2, median)
#         stab[-(i)] <- (colSums(tmp$t!=0)/tmp$R)[-1]
#         c(theta, stab)
#     }
#     list(a=res[,1], b=postProcess(res[,1:p+1]), stab=res[,(p+2):ncol(res)])
# }

##' @title func.E.stab
##'
##' @param dat.tss relative abundances matrix (each OTU in one row)
##' @param sample.filter filter out samples contain outliers in Y
##' @param m estimated biomass values
##' @param ncpu number of CPUs (default: 4)
##' @param perc percentage of samples to take for each iteration
##' @param niter number of iterations to run
##' @param ... additional parameters for beemStatic:::infer
##' @importFrom doMC registerDoMC
##' @import foreach
##' @description E-step with subsampling to estimate stability
##' @author Chenhao Li, Niranjan Nagarajan
func.E.stab <- function(dat.tss, m, sample.filter, ncpu=4, perc=0.6, niter=100, ...){
    registerDoMC(ncpu)
    p <- nrow(dat.tss)
    res <- foreach(i=1:p, .combine=rbind) %do% {
        message(paste0("Resampling for species ", i))
        fil <- dat.tss[i,]!=0 & !sample.filter[i,]
        X <- t(rbind(1/m, dat.tss[-i,])[,fil])
        Y <- dat.tss[i, fil]
        theta <- rep(0, p+1)
        stab <- rep(1, p)
        theta[i+1] <- -1 ## -beta_{ii}/beta_{ii}
        tmp <- (sapply(1:niter, function(x) sub_stat(data=cbind(Y,X), perc) ))
        theta[-(i+1)] <- apply(tmp, 1, median)
        stab[-(i)] <- (rowSums(tmp!=0)/niter)[-1]
        c(theta, stab)
    }
    list(a=res[,1], b=postProcess(res[,1:p+1]), stab=res[,(p+2):ncol(res)])
}

##' @title sub_stat
##'
##' @param data data in the format of cbind(Y, X)
##' @param perc percentage of samples to take for each iteration
##' @param ... additional parameters for beemStatic:::infer
##' @description Resampling for the inference process
##' @author Chenhao Li, Niranjan Nagarajan
sub_stat <- function(data, perc, ...) {
    n <- nrow(data)
    indices <- sample(1:n, n*perc)
    X_s <- data[indices, -1]
    Y_s <- data[indices, 1]
    fit <- infer(Y_s, X_s)
    return(fit[1:ncol(X_s)])
}

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

##' @title entropy
##'
##' @param v a vector input (treated as binary, i.e. zero/non-zero)
##' @description Compute the entropy of a binary vector
##' @author Chenhao Li, Niranjan Nagarajan
entropy <- function(v){
    p <- sum(v == 0)/length(v)
    ifelse(p==0 || p==1, 0, -(p*log2(p) + (1-p)*log2(1-p)) )
}

##' @title func.E
##'
##' @param dat.tss relative abundances matrix (each OTU in one row)
##' @param sample.filter filter out samples contain outliers in Y
##' @param lambda.inits initial lambda values
##' @param m estimated biomass values
##' @param ncpu number of CPUs (default: 4)
##' @param center center data or not
##' @importFrom doMC registerDoMC
##' @import foreach
##' @description E-part of BEEM, estimate model parameters with inferred m
##' @author Chenhao Li, Niranjan Nagarajan
func.E <- function(dat.tss, m, sample.filter, lambda.inits=NULL, ncpu=4, center=FALSE, ...){
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
        tmp <- infer(Y,X, lambda.init=lambda.inits[i], ...)
        theta[-(i+1)] <- tmp[1:p]
        e2 <- rep(NA, ncol(dat.tss))
        e2[fil] <- tmp[(p+1):(length(tmp)-1)]
        lambda <- tmp[length(tmp)]
        c(theta, e2, lambda)
    }
    ## check if there is not enough information
    uncertain <- foreach(i=1:p, .combine=rbind) %dopar% {
        fil <- dat.tss[i,]!=0 & !sample.filter[i,]
        X <- dat.tss[, fil]
        apply(X, 1, entropy)
        ## abs(rowSums(X!=0)/ncol(X) - 0.5) ## non-zero entries
    }

    list(a=res[,1], b=postProcess(res[,1:p+1]), e2=res[,(p+2):(ncol(res)-1)],
         uncertain=uncertain, lambdas=res[,ncol(res)])
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
##' @param dev dev * IQR from median will be filtered out (default: 0, nothing to remove)
##' @description pre-process data
##' @author Chenhao Li, Niranjan Nagarajan
preProcess <- function(dat, dev=0){
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
##' @importFrom reshape2 melt
##' @param a scaled growth rates
##' @param b scaled interaction matrix
##' @param vnames variable names
##' @description Function to convert parameter vector a and matrix b to MDSINE's output format
##' @author Chenhao Li, Niranjan Nagarajan
formatOutput <- function(a, b, vnames){
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
    score <- (err - median(err, na.rm = TRUE))/IQR(err, na.rm = TRUE)
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
##' @param m.init initial biomass values (default: use CSS normalization)
##' @param max.iter maximal number of iterations (default 30)
##' @param warm.iter number of iterations to run before removing any samples (default: run until convergence and start to remove samples)
##' @param resample number of iterations to resample the data to compute stability of the interaction parameters (default: 0 - no resampling)
##' @param alpha the alpha parameter for the Elastic Net model (1-LASSO [default], 0-RIDGE)
##' @param refresh.iter refresh the removed samples every X iterations (default: 3)
##' @param verbose print out messages
##' @param debug output debugging information (default FALSE)
##' @description Iteratively estimating scaled parameters and biomass
##' @export
##' @author Chenhao Li, Niranjan Nagarajan
func.EM <- function(dat, ncpu=4, scaling=1000, dev=Inf, m.init=NULL,
                    max.iter=30, warm.iter=NULL, resample=0,
                    alpha=1, refresh.iter=3, debug=FALSE, verbose=TRUE){

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
    if(is.null(m.init)) {
        m.iter <- scaling * tmp/median(tmp)
    }else{
        m.iter <- scaling * m.init/median(m.init)
    }
    trace.m <- matrix(m.iter)
    trace.p <- apply(expand.grid(spNames, spNames), 1, function(x) paste0(x[2], '->', x[1]))
    trace.p <- c(paste0(NA, '->',spNames), trace.p)
    trace.p <- data.frame(name=trace.p)
    trace.lambda <- data.frame(spNames)
    lambda.inits <- NULL
    ## flags
    remove_non_eq <- FALSE
    removeIter <- 0

    ## constants
    m1 <- matrix(rep(1,nrow(sample.filter.iter)))
    ## EM
    for(iter in 1:max.iter){
        if(verbose) message(paste0("############# Run for iteration ", iter,": #############"))
        if(verbose) message("E-step: estimating scaled parameters...")
        if(iter>=2) lambda.inits <- lambdas
        tmp.p <- func.E(dat.tss, m.iter, sample.filter.iter, ncpu=ncpu, method='glmnet',
                        lambda.inits=lambda.inits, alpha=alpha)
        err.p <- tmp.p$e2
        lambdas <- tmp.p$lambdas
        uncertain <- tmp.p$uncertain
        if(debug){
            print(rowSums(tmp.p$b!=0))
            # par(mfrow=c(6,6),mar=c(1,1,1,1))
            # for(i in 1:nrow(dat.tss)){
            #     plot(as.numeric(err.p[i,])~as.numeric(dat.tss[i,]))
            # }
            ##par(mfrow=c(1,2))
            plot(1-rowMeans(err.p, na.rm = TRUE)/apply(dat.tss[,], 1, function(x) var(x[x!=0], na.rm=TRUE) ))
            ##abline(h=0.5)
        }
        if(verbose) message("M-step: estimating biomass...")
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
            if(verbose) message(paste0("Number of samples removed (detected to be non-static): ",
                           sum(colSums(sample.filter.iter)>0)))
            removeIter <- removeIter + 1
        }

        m.iter <- m.iter*scaling/median(m.iter[colSums(sample.filter.iter)==0], na.rm=TRUE)
        trace.m <- cbind(trace.m, m.iter)
        trace.p <- cbind(trace.p, formatOutput(tmp.p$a, tmp.p$b, spNames)$value)
        trace.lambda <- cbind(trace.lambda, lambdas)
        criterion <- median(abs((trace.m[, iter+1] - trace.m[, iter])/trace.m[, iter]))<1e-3
        if (!is.null(warm.iter) && iter > warm.iter && !remove_non_eq){
            if(verbose) message("Start to detect and remove bad samples...")
            remove_non_eq <- TRUE
        }
        if (iter > 5 && !remove_non_eq && criterion && is.finite(dev)) {
            if(verbose) message("Converged and start to detect and remove bad samples...")
            remove_non_eq <- TRUE
        }
        if (((removeIter > 5 && remove_non_eq) || is.infinite(dev)) && criterion) {
            if(verbose) message("Converged!")
            break
        }
    }
    res.resample <- NULL
    if(resample!=0){
        if(verbose) message("Estimating stability...")
        ##res.resample <- func.E.stab(dat.tss, m.iter, sample.filter.iter, ncpu=ncpu, alpha=alpha)
        res.resample <- resample.EM(dat[, !colSums(sample.filter.iter) > 0], m=m.iter[!colSums(sample.filter.iter) > 0],
                                    perc=0.6, res.iter=resample,
                                    ncpu=ncpu, scaling=scaling, dev=dev, refresh.iter=refresh.iter, alpha=alpha,
                                    max.iter=20, warm.iter=0, resample=0, debug=FALSE, verbose=FALSE)
        res.resample$a.summary <- apply(res.resample$res.a, 1, median)
        res.resample$b.summary <- matrix(apply(res.resample$res.b, 1, function(x) median(x)), nrow(dat))
        res.resample$b.stab <- matrix(rowSums(res.resample$res.b != 0)/ncol(res.resample$res.a), nrow(dat))
    }
    list(trace.m=trace.m, trace.p=trace.p, err.m=err.m, err.p=err.p,
         b.uncertain = uncertain, trace.lambda=trace.lambda,
         resample=res.resample,
         sample2rm = which(colSums(sample.filter.iter) > 0))
}

##' @title resample.EM
##'
##' @param data data in the format of cbind(Y, X)
##' @param m biomass initialization
##' @param perc percentage of samples to take for each iteration
##' @param res.iter number of resample iteration
##' @param ... additional parameters for beemStatic::func.EM
##' @description Resampling for the inference process
##' @export
##' @author Chenhao Li, Niranjan Nagarajan
resample.EM <- function(data, m, perc, res.iter, ...) {
    n <- ncol(data)
    p <- nrow(data)
    res <- foreach(i=1:res.iter, .combine=cbind) %do% {
        message(paste0("#### Resample iteration: ", i, " #####"))
        indices <- sort(sample(1:n, n*perc))
        tmp <- func.EM(data[, indices], m.init=m[indices], ...)
        tmp$trace.p[, ncol(tmp$trace.p)]
    }
    list(res.a = res[1:p,], res.b = res[-c(1:p),])
}

