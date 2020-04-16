##' @title inference
##'
##' @importFrom selectiveInference fs
##' @importFrom selectiveInference fsInf
##' @param dat original data used for BEEM-static prediction
##' @param beem BEEM-static ouput
##' @param ncpu number of CPUs (default: 1)
##' @param p.values report p-values computed from the confidence scores (default: FALSE)
##' @description Infer confidence based on stepwise variable selection
##' @export
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
inference <- function(dat, beem,  ncpu=1, p.values=FALSE){
  ## currently doesn't support perturbation parameters
  registerDoParallel(ncpu)
  m <- beem2biomass(beem)
  dat.tss <- tss(dat)
  if(length(beem$sample2rm)>0){
    dat.tss <- dat.tss[, -beem$sample2rm]
    m <- m[-beem$sample2rm]
  }
  p <- nrow(dat.tss)
  param <- beem2param(beem)
  res <- foreach(i=1:p, .combine=rbind) %dopar% {
    fil <- dat.tss[i,]!=0 #& !sample.filter[i,]
    X <- t(rbind(1/m, dat.tss[-i,])[,fil])
    Y <- dat.tss[i, fil]
    tvs <- rep(0, p)
    pvs <- rep(1,p)
    tvs[i] <- NA ## -beta_{ii}/beta_{ii}
    pvs[i] <- NA
    beta <- c(param$b.est[i,-i])
    if (any(beta!=0)){
      id.non0 <- colSums(X!=0) > 1
      fit <- fs(X[, id.non0], Y, intercept = F)
      x <- fsInf(fit,type="aic")
      z <- abs(x$sign * x$vmat %*% x$y/(x$sigma * sqrt(rowSums(x$vmat^2))))[,1]
      ids <- x$vars - 1
      pvs[-i][id.non0[-1]][ids[ids!=0]] <- x$pv[ids!=0]
      tvs[-i][id.non0[-1]][ids[ids!=0]] <- z[ids!=0]
    }
    c(tvs, pvs)
  }
  if (!p.values) {
      return(res[, 1:p])
  }
  return(list(conf=res[,1:p], p.value=res[,(p+1):(2*p)]))
}

##' @title predict.beem
##'
##' @param dat.new OTU count/relative abundance matrix (each OTU in one row)
##' @param object output of the EM algorithm
##' @param dev deviation of the error (for one sample) from the model to be excluded
##' @param ncpu number of CPUs (default: 4)
##' @param pert.new external perturbation presence matrix (each perturbation in one row, each sample in one column) (Default: NULL)
##' @importFrom doParallel registerDoParallel
##' @description Use a trained BEEM-static model to predict biomass, deviation from steady states and violation of model assumption
##' @export
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
predict.beem <- function(object, dat.new, dev, ncpu=4, pert.new = NULL){
  ### currently not ready for an S3method yet
  beem <- object
  param <- beem2param(beem)
  dat.new.tss <- tss(dat.new)
  if (!is.null(pert.new)){
    tmp.m <- func.M(dat.new.tss, param$a.est, param$b.est, c = param$c.est, perturbation.presence = pert.new, ncpu = 1)
  } else {
    tmp.m <- func.M(dat.new.tss, param$a.est, param$b.est, ncpu=ncpu)
  }
  m.pred <- tmp.m[,1]
  err.m.pred <- tmp.m[,-1]

  registerDoParallel(ncpu)
  p <- nrow(dat.new.tss)
  e <- foreach(i=1:p, .combine=rbind) %dopar% {
    X <- t(rbind(1/m.pred, dat.new.tss[-i,]))
    Y <- dat.new.tss[i, ]
    coefs <- c(param$a.est[i], param$b.est[i,-i])
    if(!is.null(param$c.est)){
      pert <- t(pert.new)
      coefs_pert <- param$c.est[i,]
      as.numeric(Y - (X %*% coefs)[,1] - (1/m.pred) * ((pert %*% coefs_pert)[,1]) )
    }else{
      as.numeric((Y - (X %*% coefs)[,1] ) )
    }
  }
  ## predicting the equilibrium based on the current results
  eq.pred <- foreach(i=1:ncol(dat.new), .combine = cbind) %dopar%{
    sel <- dat.new[,i]>0 ## species present
    tmp <- matrix(0, length(sel), 1)
    tmp[sel] <- - solve(param$b.est[sel,sel]) %*% (param$a.est/m.pred[i])[sel]
    tmp
  }
  tmp.err <- apply(abs(beem$err.p)/tss(beem$dat), 2, median, na.rm=TRUE)
  score <- (abs(e)/dat.new.tss - median(tmp.err, na.rm = TRUE))/IQR(tmp.err, na.rm = TRUE)
  score[is.na(score)] <- Inf
  isBad <- score > dev
  return(list(biomass.pred=m.pred, dev.from.eq=t(err.m.pred), isBad=isBad, eq.pred=eq.pred))
}

##' @title print.beem
##'
##' @param object output of the EM algorithm
##' @description print method for beem object
##' @export
##' @author Chenhao Li
print.beem <- function(object){
    cat("\nEM-call: \n")
    print(object$call)
    cat("\n\n")
    invisible(object)
}
