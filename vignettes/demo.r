############################################################################
## Preprocessing
#### Generate simulated data with noise
df <- read.table('../data/demo.n500_p20_eq1.counts.txt', head=F, row.names=1, sep='\t')

aux.poisson <- function(v, depth){
    v <- round(v/sum(v) * depth)
    sapply(v, function(x) rpois(1, x))
}
df.noise <- apply(df, 2, aux.poisson, depth=5000)

params <- read.table('../data/demo.n500_p20_eq1.truth.txt', head=T)

truth2param <- function(truth){
    a.truth <- truth[is.na(truth$source_taxon),3]
    p <- length(a.truth)
    b.truth <- matrix(truth[!is.na(truth$source_taxon),3],p,p)
    tmp.diag <- truth[truth$source_taxon==truth$target_taxon,3][-(1:p)]
    a.truth <- -a.truth/tmp.diag
    b.truth <- -b.truth/tmp.diag
    return(list(a.truth=a.truth, b.truth=b.truth))
}

auc.b <- function(b.est, b.true){
    diag(b.est) <- NA
    diag(b.true) <- NA
    est <- c(b.est)
    lab <- c(b.true)
    est <- abs(est[!is.na(est)])
    lab <- lab[!is.na(lab)]
    lab <- (lab!=0 ) *1
    pROC::plot.roc(lab,est, print.auc = TRUE)
}

scaled.params <- truth2param(params)

out <- list(dat.w.noise=df.noise, dat.wo.noise=df,
            scaled.params=scaled.params,
            original.params=params, biomass.true=colSums(df))

saveRDS(out, '../data/demo_dat.rds')

##############################################################################

dat <- readRDS('../data/demo_dat.rds')
attach(dat)

res <- func.EM(dat.w.noise, ncpu=10, scaling=median(dat$biomass.true), max.iter=40)

diagnoseBiomass(res, true.biomass = biomass.true)

plot(res$trace.m[,30], biomass.true, xlab='BEEM biomass estimation', ylab='True biomass')

est <- beem2param(res)

plot(est$a.est, scaled.params$a.truth)

plot(est$b.est, scaled.params$b.truth)


library(SpiecEasi)

se <- spiec.easi(t(dat.w.noise), method='mb')
graph <- as.matrix(getOptMerge(se))

auc.b(graph, scaled.params$b.truth)
auc.b(est$b.est, scaled.params$b.truth)
