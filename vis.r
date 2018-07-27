##' @title diagnoseBiomass
##' 
##' @param beem.out output of a beem run
##' @param true.biomass measured/true biomass (default: the biomass of last iteration)
##' @param alpha transparency parameter of the biomass lines
##' @description plot the trace of the fitted error of biomass
##' @author Chenhao Li, Niranjan Nagarajan
##' @export
diagnoseBiomass <- function(beem.out, true.biomass=NA, alpha=0.1,...){
    trace.m <- beem.out$trace.m
    nIter <- ncol(trace.m)
    if(any(is.na(true.biomass))){
        true.m <- trace.m[,nIter]
    }else{
        true.m <- true.biomass
    }
    rel.err.m <- (t((trace.m-true.m)/true.m))*100
    col <- rep(rgb(0,0,0,alpha), nrow(trace.m))
    col[beem.out$sample2rm] <- rgb(1,1,0, alpha)
    matplot(rel.err.m, type='l',
            xlab="Iterations",
            ylab="Relative difference (%)",
            main='Biomass trace',
            lty=1, lwd = 3,
            col=col,
            ...
            )
    lines(x=1:ncol(trace.m),y=apply(rel.err.m,1,median), col='red', lwd=5)    
}

##' @title pcoa
##' 
##' @param countData OTU/species abundance table (each row is one species, each column is one site)
##' @description perform a PCoA analysis using Bray-Curtis distance
##' @importFrom vegan vegdist
##' @import ggplot2
##' @author Chenhao Li, Niranjan Nagarajan
##' @export
pcoa <- function(countData){
    dat.pcoa <- cmdscale(vegan::vegdist(t(countData)), eig=TRUE)
    per.var <- (dat.pcoa$eig/sum(dat.pcoa$eig))[1:2] * 100
    ggplot(data.frame(dat.pcoa$points), aes(x=X1, y=X2)) +
        geom_point(size=2) +
        labs(x=paste0('CMD1 (' ,round(per.var[1], 2),'%)'),
             y=paste0('CMD2 (',round(per.var[2], 2),'%)'))    
}

