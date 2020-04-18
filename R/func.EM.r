##' @title detectBadSamples
##'
##' @param err the errors estimated from regression
##' @param threshold threshold to filter out samples
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
detectBadSamples <- function(err, threshold){
    score <- (err - median(err, na.rm = TRUE))/IQR(err, na.rm = TRUE)
    score[is.na(score)] <- Inf
    return(score > threshold)
}


##' @title func.EM
##'
##' @param dat OTU count/relative abundance matrix (each OTU in one row)
##' @param external.perturbation external perturbation presence matrix (each perturbation in one row, each sample in one column) (Default: NULL)
##' @param ncpu number of CPUs (default: 1)
##' @param scaling a scaling factor to keep the median of all biomass constant (default: 1000)
##' @param dev deviation of the error (for one sample) from the model to be excluded (default: Inf - all the samples will be considered)
##' @param m.init initial biomass values (default: use CSS normalization)
##' @param max.iter maximal number of iterations (default 30)
##' @param lambda.iter number of iterations to run before fixing lambda (default: 2)
##' @param warm.iter number of iterations to run before removing any samples (default: run until convergence and start to remove samples)
##' @param lambda.choice 1: use lambda.1se for LASSO, 2: use lambda.min for LASSO, a number between (0, 1): this will select a lambda according to (1-lambda.choice)*lambda.min + lambda.choice*lambda.1se
##' @param resample number of iterations to resample the data to compute stability of the interaction parameters (default: 0 - no resampling)
##' @param alpha the alpha parameter for the Elastic Net model (1-LASSO [default], 0-RIDGE)
##' @param refresh.iter refresh the removed samples every X iterations (default: 1)
##' @param epsilon convergence threshold (in relative difference): uqn of the relative error in biomass changes (default 1e-3)
##' @param verbose print out messages
##' @param debug output debugging information (default FALSE)
##' @description Iteratively estimating scaled parameters and biomass
##' @export
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
func.EM <- function(dat, external.perturbation = NULL, ncpu=1, scaling=1000, dev=Inf, m.init=NULL,
                    max.iter=30, lambda.iter=2, warm.iter=NULL, lambda.choice=1, resample=0,
                    alpha=1, refresh.iter=1, epsilon=1e-3,
                    debug=FALSE, verbose=TRUE){

    cl <- match.call()
    ## pre-processing
    dat.init <- preProcess(dat, dev=0)
    dat.tss <- dat.init$tss
    spNames <- rownames(dat)

    ## ensure valid samples
    temp <- colSums(dat.tss, na.rm = TRUE) == 0
    if(any(temp)){
        stop(paste0('Sample ', which(temp), ' has zero total abudance...\n'))
    }
    ## ensure valid variables
    temp <- rowSums(dat.tss>0) <= 10
    if(any(temp)){
        stop(paste0('Variable ', which(temp), ' has too many low abundance values...\n'))
    }
    ## ensure parameters are put correctly
    if(lambda.iter < 2){
        stop('The value of lambda.iter should not be smaller than 2...')
    }
    ## initialization
    sample.filter.iter <- dat.init$sample.filter
    tmp <- css(t(dat.tss))$normFactors
    #tmp <- biomassInit(dat.tss)
    if(is.null(m.init)) {
        m.iter <- scaling * tmp/median(tmp)
    }else{
        m.iter <- scaling * m.init/median(m.init)
    }
    trace.m <- matrix(m.iter)
    trace.p <- apply(expand.grid(spNames, spNames), 1, function(x) paste0(x[2], '->', x[1]))
    trace.p <- c(paste0(NA, '->',spNames), trace.p)
    if (!is.null(external.perturbation)) {
        external.perturbation.Names <- rownames(external.perturbation)
        trace.p.external.perturbation <- apply(expand.grid(spNames, external.perturbation.Names), 1, function(x) paste0(x[2], '->', x[1]))
        trace.p <- append(trace.p, trace.p.external.perturbation)
    }
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
        if(iter>=lambda.iter) lambda.inits <- lambdas
        if(is.null(external.perturbation)){
            ext.pert <- NULL
        }else{
            ext.pert <- 1/m.iter * external.perturbation
        }
        tmp.p <- func.E(dat.tss, external.perturbation = ext.pert,
                        m.iter, sample.filter.iter, ncpu=ncpu, method='glmnet',
                        lambda.inits=lambda.inits, alpha=alpha)
        err.p <- tmp.p$e
        lambdas <- tmp.p$lambdas
        uncertain <- tmp.p$uncertain
        if(debug){
            print(rowSums(tmp.p$b!=0))
            plot(1-rowMeans(err.p^2, na.rm = TRUE)/apply(dat.tss[,], 1, function(x) var(x[x!=0], na.rm=TRUE) ))
            ##abline(h=0.5)
        }

        if(verbose) message("M-step: estimating biomass...")
        if (!is.null(external.perturbation)) {
            tmp.m <- func.M(dat.tss, tmp.p$a, tmp.p$b, c = tmp.p$perturbation.coefficients, perturbation.presence = external.perturbation, ncpu)
        } else {
            tmp.m <- func.M(dat.tss, tmp.p$a, tmp.p$b, c = NULL, perturbation.presence = NULL, ncpu)
        }
        m.iter <- tmp.m[,1]
        err.m <- tmp.m[,-1]

        if(remove_non_eq){
            ## clear up removed samples every X iterations
            if (iter %% refresh.iter == 0 ) {
                sample.filter.iter <- dat.init$sample.filter
            }
            bad.samples <- detectBadSamples(apply(abs(t(err.p/dat.tss))/m.iter, 1, median, na.rm=TRUE), dev)
            #plot(apply(abs(t(err.p/dat.tss))/m.iter, 1, median, na.rm=TRUE), col=ifelse(bad.samples, 'red', 'black'))
            #bad.samples <- detectBadSamples(apply(abs(err.p)/dat.tss, 2, median, na.rm = TRUE), dev)
            sample.filter.iter <- (m1 %*% bad.samples)>0  | sample.filter.iter
            if(verbose) message(paste0("Number of samples removed (detected to be non-static): ",
                                       sum(colSums(sample.filter.iter)>0)))
            removeIter <- removeIter + 1
        }
        m.fil <- colSums(sample.filter.iter)==0
        m.iter <- m.iter*scaling/median(m.iter, na.rm=TRUE)
        trace.m <- cbind(trace.m, m.iter)
        if (!is.null(external.perturbation)) {
            trace.p <- cbind(trace.p, formatOutput(tmp.p$a, tmp.p$b, c = tmp.p$perturbation.coefficients,
                                                   spNames, pnames = external.perturbation.Names)$value)
        } else {
            trace.p <- cbind(trace.p, formatOutput(tmp.p$a, tmp.p$b, c = NULL, spNames, pnames = NULL)$value)
        }

        trace.lambda <- cbind(trace.lambda, lambdas)
        summary.m <- apply(trace.m[m.fil, ], 1, function(x) median(rev(x)[1:min(5,length(x))]))
        #criterion <- quantile(abs((trace.m[m.fil, iter+1] - summary.m)/summary.m), 1)<epsilon
        criterion <- median(abs((trace.m[m.fil, iter+1] - trace.m[m.fil, iter])/trace.m[m.fil, iter]))<epsilon #0.0025

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
        if (!is.null(external.perturbation)) {
            res.resample <- resample.EM(dat[, m.fil], external.perturbation = external.perturbation[, m.fil], m=m.iter[m.fil],
                                        perc=0.6, res.iter=resample,
                                        ncpu=ncpu, scaling=scaling, dev=dev, refresh.iter=refresh.iter, alpha=alpha,
                                        max.iter=20, warm.iter=0, resample=0, debug=FALSE, verbose=FALSE)
        } else {
            res.resample <- resample.EM(dat[, m.fil], external.perturbation = NULL, m=m.iter[m.fil],
                                        perc=0.6, res.iter=resample,
                                        ncpu=ncpu, scaling=scaling, dev=dev, refresh.iter=refresh.iter, alpha=alpha,
                                        max.iter=20, warm.iter=0, resample=0, debug=FALSE, verbose=FALSE)
        }
        res.resample$a.summary <- apply(res.resample$res.a, 1, median)
        res.resample$b.summary <- matrix(apply(res.resample$res.b, 1, function(x) median(x)), nrow(dat))
        res.resample$b.stab <- matrix(rowSums(res.resample$res.b != 0)/ncol(res.resample$res.a), nrow(dat))
        if (!is.null(external.perturbation)) {
            res.resample$c.summary <- apply(res.resample$res.c, 1, function(x) median(x))
            res.resample$c.stab <- matrix(rowSums(res.resample$res.c != 0)/ncol(res.resample$res.a), nrow(dat))
        }
    }
    structure(
        list(trace.m=trace.m, trace.p=trace.p, trace.lambda=trace.lambda,
             err.m=err.m, err.p=err.p,
             b.uncertain = uncertain,
             #resample=res.resample,
             dev=dev,
             sample2rm = which(!m.fil),
             input=dat,
             call=cl),
        class="beem"
    )
}



##' @title formatOutput
##' @importFrom reshape2 melt
##' @param a scaled growth rates
##' @param b scaled interaction matrix
##' @param c scaled external perturbation matrix (default: NULL)
##' @param vnames species names
##' @param pnames external perturbation names (default: NULL)
##' @description Function to convert parameter vector a and matrix b to MDSINE's output format
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
formatOutput <- function(a, b, c = NULL, vnames, pnames = NULL){
    p <- length(a)
    param <- data.frame(parameter_type=c(rep("growth_rate",p), rep("interaction",p*p)))
    param$source_taxon <- c(rep(NA,p),rep(vnames, each = p))
    param$target_taxon <- vnames
    param$value <- 0
    param$significance <- NA
    tmp.b <- melt(t(b))
    tmp.b <- tmp.b[order(tmp.b$Var1),]
    param$value <- c(a, tmp.b$value)
    if (!is.null(c) || !is.null(pnames)) {
        q <- ncol(c) #Finding the number of external perturbationss
        ext.perturb <- data.frame(parameter_type=c(rep("external perturbations", p*q)))
        ext.perturb$source_taxon <- c(rep(NA, p*q))
        ext.perturb$target_taxon <- vnames
        ext.perturb$value <- 0
        ext.perturb$significance <- NA
        tmp.c <- melt(c)
        ext.perturb$value <- c(tmp.c$value)
        param <- rbind(param, ext.perturb)
    }
    param$external_perturbation <- c(rep(NA,p), rep(NA,p*p), rep(pnames, each = p)) #Creating the vector for external perturbations
    param
}
