##' @title beem2param
##'
##' @param beem a BEEM object
##' @description extract parameter estimates from a BEEM object
##' @export
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
beem2param <- function(beem){
    p <- nrow(beem$input) #Number of species
    num.perturb <- (nrow(beem$trace.p) - p - p*p)/p
    tmp <- beem$trace.p[, ncol(beem$trace.p)]
    a.est <- tmp[1:p]
    b.est <- matrix(tmp[-c(1:p)], p, p)
    if (num.perturb == 0) {
        return (list(a.est=a.est, b.est=b.est))
    } else {
        c.est <- matrix(tmp[(nrow(beem$trace.p) - (num.perturb*p - 1)):nrow(beem$trace.p)],
                        ncol = num.perturb)
        return (list(a.est=a.est, b.est=b.est, c.est=c.est))
    }
}

##' @title beem2biomass
##'
##' @param beem a BEEM object
##' @description extract biomass estimates from a BEEM object
##' @export
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
beem2biomass <- function(beem){
    niter <- ncol(beem$trace.m)
    return (beem$trace.m[,niter])
}
