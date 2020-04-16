##' @title func.E
##'
##'
##' @param dat.tss relative abundances matrix (each OTU in one row)
##' @param external.perturbation external perturbation presence matrix *1/m (each perturbation in one row, each sample in one column) (Default: NULL)
##' @param sample.filter filter out samples contain outliers in Y
##' @param lambda.inits initial lambda values
##' @param m estimated biomass values (1 X no. of samples) matrix
##' @param ncpu number of CPUs (default: 4)
##' @param center center data or not
##' @param ... additional parameters for `beemStatic::infer`
##' @importFrom doParallel registerDoParallel
##' @import foreach
##' @description E-part of BEEM, estimate model parameters with inferred m
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
func.E <- function(dat.tss, external.perturbation = NULL, m, sample.filter, lambda.inits=NULL, ncpu=4, center=FALSE, ...){
    ## infer parameter for each OTU
    registerDoParallel(ncpu)
    p <- nrow(dat.tss)
    if (!is.null(external.perturbation)) {
        k <- nrow(external.perturbation) #Number of external perturbations
    }
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
        if (!is.null(external.perturbation)) {
            perturbation.coefficients <- rep(0, k) #Creating a vector to store perturbation coefficients
            external.perturbation.filtered <- external.perturbation[,fil] #Filtering out the external perturbations corresponding to samples which were filtered out
            X <- cbind(X, t(external.perturbation.filtered))
        }

        tmp <- infer(Y,X, lambda.init=lambda.inits[i], ...)

        if (!is.null(external.perturbation)) {
            theta[-(i+1)] <- tmp[1:p]
            perturbation.coefficient <- tmp[seq(p+1,p+k)]
            e <- rep(NA, ncol(dat.tss))
            e[fil] <- tmp[(p+k+1):(length(tmp)-1)]
            lambda <- tmp[length(tmp)]
            c(theta, perturbation.coefficient, e, lambda)
        } else {
            theta[-(i+1)] <- tmp[1:p]
            e <- rep(NA, ncol(dat.tss))
            e[fil] <- tmp[(p+1):(length(tmp)-1)]
            lambda <- tmp[length(tmp)]
            c(theta, e, lambda)
        }
    }

    ## check if there is not enough information
    uncertain <- foreach(i=1:p, .combine=rbind) %dopar% {
        fil <- dat.tss[i,]!=0
        X <- dat.tss[, fil]
        apply(X, 1, entropy)
        ## abs(rowSums(X!=0)/ncol(X) - 0.5) ## non-zero entries
    }
    if (!is.null(external.perturbation)) {
        list(a=res[,1], b=postProcess_species(res[,1:p+1]),
             perturbation.coefficients = postProcess_perturbation(matrix(res[,seq(p+2,p+k+1)], ncol = k), res[,1]),
             e=res[,(p+k+2):(ncol(res)-1)],
             uncertain=uncertain, lambdas=res[,ncol(res)])
    } else {
        list(a=res[,1], b=postProcess_species(res[,1:p+1]), e=res[,(p+2):(ncol(res)-1)],
             uncertain=uncertain, lambdas=res[,ncol(res)])
    }
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
##' @param lambda.adjust.up allow adjusting lambda (increase) by 50\% in each iteration
##' @param lambda.adjust.down allow adjusting lambda (decrease) by 50\% in each iteration
##' @param nfolds number of folds for glmnet cross-validation
##' @import glmnet
##' @description Infer parameters with blasso
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
infer <- function(Y, X, method='glmnet', intercept=FALSE, seed=0, alpha=1, nfolds=3,
                  lambda.init=NULL, lambda.choice=1,
                  lambda.adjust.up=TRUE, lambda.adjust.down=TRUE){
    lambda.lower <- 1e-9
    lambda.upper <- 1
    set.seed(seed)
    if(method=='lm'){
        return(as.numeric(coef(lm(Y~X+0))))
    }
    if(method=='glmnet'){
        penalty <- c(0,rep(1, ncol(X)-1))
        maxmin <- median(Y) + 5 * IQR(Y) %*% c(1,-1)
        idx <- Y <= maxmin[1] & Y >= maxmin[2]
        if(is.null(lambda.init)){
            ## select a lambda range first
            lambda.init <- rev(exp(seq(log(lambda.lower), log(lambda.upper), length.out = 50)))
            fit <- cv.glmnet(X[idx,], Y[idx], intercept=intercept, lambda=lambda.init, nfolds=nfolds,
                             penalty.factor=penalty, alpha=alpha)
            lambda <- rev(exp(seq(
                max(log(fit$lambda.min/10), log(lambda.lower)),
                min(log(fit$lambda.min*10), log(lambda.upper)),
                length.out = 50)))
        }else{
            ## adjust lambda by up to 50%
            lambda <- rev(exp(seq(
                max(log(lambda.init/2), log(lambda.lower)), ## bounded by the min lambda
                min(log(lambda.init*1.5), log(lambda.upper)), ## bounded by the max lambda
                length.out = 20)))
        }

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
        e <- as.numeric((Y-(X %*% coefs)[,1]))
        rm(.Random.seed, envir=.GlobalEnv)
        return(c(coefs, e, s))
    }
}

##' @title postProcess_species
##'
##' @param b scaled interaction matrix (diagonal/self-interaction is -1)
##' @param dev scale of the interaction < dev * 1, the value will be set to 0
##' @description post-process parameters to remove small entries
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
postProcess_species <- function(b, dev=1e-5){
    b[abs(b)<dev] <- 0
    b
}

##' @title postProcess_perturbation
##'
##' @param b scaled perturbation matrix (comparing to growthrate)
##' @param dev scale of the interaction < dev * 1, the value will be set to 0
##' @description post-process parameters to remove small entries
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
postProcess_perturbation <- function(c, growth.rate, dev=1e-5){
    c[abs(c/growth.rate)<dev] <- 0
    c
}

##' @title entropy
##'
##' @param v a vector input (treated as binary, i.e. zero/non-zero)
##' @description Compute the entropy of a binary vector
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
entropy <- function(v){
    p <- sum(v == 0)/length(v)
    ifelse(p==0 || p==1, 0, -(p*log2(p) + (1-p)*log2(1-p)) )
}
