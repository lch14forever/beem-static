##' @title preProcess
##'
##' @param dat OTU count/relative abundance matrix (each OTU in one row)
##' @param dev dev * IQR from median will be filtered out (default: 0, nothing to remove)
##' @description pre-process data
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
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

##' @title tss
##'
##' @param x count matrix (each OTU in one row)
##' @description Calculate the relative abundances by dividing the total abundance (Total Sum Sacling)
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
tss <- function(x){ apply(x, 2, function(x) x/sum(x)) }


##' @title css
##'
##' @param m matrix of data (variables in columns, measurements in rows)
##' @param p quantile used for normalization (default: 0.5)
##' @description Function to perform cumulative sum scaling (CSS)
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
css <- function(m, p=0.5){
    m[m==0] <- NA
    ## find the quantile in each sample
    quant <- apply(m,1,function(x) quantile(x, p=p, na.rm = TRUE))
    ## calculate normalization factor
    f <- rowSums(m*sweep(m, 1, quant, '<='),na.rm = TRUE)
    nf <- f/exp(mean(log(f)))
    dat.css <- sweep(m,1,nf, '/')
    return(list("normCounts" = dat.css, "normFactors"=1/nf))
}

##' @title biomassInit
##'
##' @param m matrix of data (variables in columns, measurements in rows)
##' @description Function to initialize biomass. This normalization assumes that all samples are at the equilibrium and there is no interaction at all.
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
biomassInit <- function(m){
    rowMed <- apply(tss(m), 1, function(x) median(x[x!=0]))
    biomass <- colSums((m > 0)*rowMed)
    return(biomass)
}
