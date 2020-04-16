##' @title func.M
##'
##' @param dat.tss relative abundance matrix (each OTU in one row)
##' @param a estimated growth rate values (scaled by self-interaction)
##' @param b estimated interaction matrix (scaled by self-interaction)
##' @param c estimated external perturbation effect matrix (scaled by self-interaction)
##' @param perturbation.presence estimated external perturbation presence (not to be multiplied by 1/m) (each perturbation in one row, each sample in one column)
##' @param ncpu number of CPUs (default: 4)
##' @importFrom doParallel registerDoParallel
##' @import foreach
##' @description M-part of BEEM, estimate biomass with inferred model parameters
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
func.M <- function(dat.tss, a, b, c = NULL, perturbation.presence=NULL, ncpu=4){
    registerDoParallel(ncpu)
    foreach(i=1:ncol(dat.tss), .combine=rbind) %dopar%{
        x <- dat.tss[,i]
        if (!is.null(c) || !is.null(perturbation.presence)) {
            s <- perturbation.presence[,i]
            norm(a, b, c = c, s = s, x)
        } else {
            norm(a, b, c = NULL, s = NULL, x)
        }
    }
}

##' @title norm
##'
##' @param a estimated growth rate values (scaled by self-interaction)
##' @param b estimated interaction matrix (scaled by self-interaction)
##' @param c estimated external perturbation effect matrix (scaled by self-interaction) (default: NULL)
##' @param s external perturbation presence vector in one sample
##' @param x relative abundances in one sample
##' @description estimate biomass with linear regression
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
norm <- function(a, b, c = NULL, s = NULL, x){
    if (!is.null(c) || !is.null(s)) {
        perturbation.eff <- c %*% s
        res <- -(a + perturbation.eff) / (b %*% x)
    } else {
        res <- -a / (b %*% x)
    }
    res <- res[x!=0]
    if(all(res < 0)){
        m <- abs(max(res))
    } else{
        m <- median(res[res>0])
    }
    if (!is.null(c) || !is.null(s)) {
        err <- (m * (b %*% x) + a + perturbation.eff)[,1]/a
    } else {
        err <- (m * (b %*% x) + a)[,1]/a
    }
    err[x==0] <- 0
    c(m, err)

}
