############################################################################
## Preprocessing
#### Generate simulated data with noise

devtools::use_data_raw()
df <- read.table('data/demo.n500_p20_eq1.counts.txt', head=F, row.names=1, sep='\t')

aux.poisson <- function(v, depth){
    v <- round(v/sum(v) * depth)
    sapply(v, function(x) rpois(1, x))
}
df.noise <- apply(df, 2, aux.poisson, depth=5000)

params <- read.table('data/demo.n500_p20_eq1.truth.txt', head=T)

truth2param <- function(truth){
    a.truth <- truth[is.na(truth$source_taxon),3]
    p <- length(a.truth)
    b.truth <- matrix(truth[!is.na(truth$source_taxon),3],p,p)
    tmp.diag <- truth[truth$source_taxon==truth$target_taxon,3][-(1:p)]
    a.truth <- -a.truth/tmp.diag
    b.truth <- -b.truth/tmp.diag
    return(list(a.truth=a.truth, b.truth=b.truth))
}

scaled.params <- truth2param(params)

beemDemo <- list(dat.w.noise=df.noise, dat.wo.noise=df,
            scaled.params=scaled.params,
            original.params=params, biomass.true=colSums(df))

devtools::use_data(beemDemo)
devtools::load_all()

##############################################################################
library(beemStatic)
data("beemDemo")
attach(beemDemo)


res <- func.EM(dat.w.noise, ncpu=4, scaling=median(biomass.true), max.iter=200, epsilon = 1e-4, dev = 5)

diagnoseFit(res, dat.w.noise, annotate = FALSE)

diagnoseBiomass(res, true.biomass = biomass.true)

plot(beem2biomass(res), biomass.true, xlab='BEEM biomass estimation', ylab='True biomass')

est <- beem2param(res)

showInteraction(res, dat.w.noise)

par(mfrow=c(1,2))
plot(est$a.est, scaled.params$a.truth, xlab='BEEM estimation', ylab='Truth', main='Growth rates')
auc.b(est$b.est, scaled.params$b.truth, main='Interaction matrix')

auc.b(inference(dat.w.noise, res), scaled.params$b.truth)

spearman <- cor(t(dat.w.noise), method='spearman')


library(SpiecEasi)
se <- spiec.easi(t(dat.w.noise), method='mb')
se.stab <- as.matrix(getOptMerge(se))


par(mfrow=c(1,3), cex=1)
auc.b(spearman, scaled.params$b.truth, is.association = TRUE, main='Spearman correlation')
auc.b(se.stab, scaled.params$b.truth, is.association = TRUE, main='SPIEC-EASI')
auc.b(est$b.est, scaled.params$b.truth, main='BEEM-static')

par(mfrow=c(1,2), cex=1)
tmp <- scaled.params$b.truth
tmp[tmp!=0 & res$b.uncertain<0.5] <- 0
auc.b(est$b.est, scaled.params$b.truth, main='BEEM-static')
auc.b(est$b.est, tmp, main='BEEM-static (excluding uncertain edges)')


