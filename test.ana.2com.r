source('/mnt/projects/lich/dream_challenge/BEEM_STATIC/dev/em_functions.r')
source('/mnt/projects/lich/dream_challenge/BEEM_STATIC/dev/evalFunc.r')

#### ad-hoc analysis
library(dplyr)

base <- '/mnt/projects/lich/dream_challenge/BEEM_STATIC/simulation/sim_p20_n400_eq/sim0103.'

df.noNoise <- read.table(paste0(base,'n400_d20.counts.txt'), head=F, row.names=1, sep='\t')

dist.to.eq <- read.table(paste0(base,'n400_d20.err_to_equil.txt'), head=F, sep='\t')
dist.to.eq[is.na(dist.to.eq)] <- max(dist.to.eq, na.rm = TRUE) * 2

matplot(abs(t(dist.to.eq)), type='b')

## library(pheatmap)
## rownames(dist.to.eq) <- colnames(dist.to.eq) <- NULL
## pheatmap(log(abs(as.matrix(dist.to.eq[,401:450]))+1), cluster_rows = F, cluster_cols = F)

truth <- read.table(paste0(base,'n400_d20.truth.txt'),head=T)


addNoise <- function(x){
    sapply(x, function(x) abs(rnorm(1, x, x*0.05)))
}

df <- apply(df.noNoise, 2, addNoise)


a.truth <- truth[is.na(truth$source_taxon),3]
p <- length(a.truth)
b.truth <- matrix(truth[!is.na(truth$source_taxon),3],p,p)
tmp.diag <- truth[truth$source_taxon==truth$target_taxon,3][-(1:p)]

a.truth <- -a.truth/tmp.diag
b.truth <- -b.truth/tmp.diag

pheatmap(b.truth, cluster_rows = F, cluster_cols = F)

ans <- formatOutput(a.truth, b.truth, vnames=rownames(df))

plot_d <- function(trace.m, trace.p){
    col = c(rep(rgb(0,0,0,0.1),400), rep(rgb(1,1,0,0.3),100) )
    par(cex.lab=1.3, mfrow = c(2,2))
    rel.err.m <- (t((trace.m-m.true)/m.true))*100
    matplot(rel.err.m, type='l',
            xlab="Iterations",
            ylab="Relative error (%)",
            main='Relative error trace',
            col=col, lty=1, lwd = 3,
            ylim=c(-100, 100)
            )
    lines(x=1:ncol(trace.m),y=apply(rel.err.m,1,median), col='red', lwd=5)
    plot(1:(ncol(trace.m)), cor(m.true, trace.m, method='spearman'), typ='b', lwd=3, pch=19,
         xlab="Iterations",
         ylab="Correlation",
         main='Correlation trace',
         )
    lines(1:ncol(trace.m), cor(m.true, trace.m), typ='b',col='blue', lwd=3, pch=19)

    plot(a.truth, trace.p[1:p, ncol(trace.p)], main='Final parameters (a)',
         xlab='Truth',ylab='BEEM', pch=19, col='red',cex=2
         )
    abline(0,1, lwd=2,lty=2)
    
    diag(b.truth) <- NA
    est <- matrix(trace.p[-c(1:p),ncol(trace.p)], p,p)
    diag(est) <- NA
    
    plot(c(b.truth),c(est) , main='Final parameters (b)',
         xlab='Truth',ylab='BEEM', pch=19, col='blue',cex=2)
    abline(0,1, lwd=2,lty=2)
}

auc.ad <- function(est, lab, ...){
    est <- abs(est[!is.na(est)])
    lab <- lab[!is.na(lab)]
    lab <- (lab!=0 ) *1
    plot.roc(lab,est, print.auc = TRUE,...)
    pROC::auc(lab,est)
}


pr.ad <- function(esti, lab, add=FALSE, ...){
    ## lab=c(b.truth)
    ##esti = c(est.regress)
    esti <- abs(esti[!is.na(esti)])

    lab <- lab[!is.na(lab)]
    lab <- (lab!=0 ) *1
    ordered <- order(esti, decreasing = TRUE)
    ## wrong signs are not considered as tps
    m <- cbind(esti, lab)[ordered,]
    n_pos <- sum((m[,2]!=0) * 1)
    n_tp <- cumsum(m[,2])
    recall <- n_tp/n_pos
    precision <- n_tp/(1:NROW(m))
    f1 <- 2*(recall*precision)/(precision+recall)
    tmp <- cbind(thre=esti[ordered], recall, precision, f1)
    if(add==FALSE){
        plot(tmp[,2], tmp[,3], ylab='Precision',xlab='Recall', type='l',
             lwd=3, cex.axis=1.6, cex.lab=1.6,...)
    }else{
        lines(tmp[,2], tmp[,3], type='l',
             lwd=3, cex.axis=1.6, cex.lab=1.6,...)
    }
    lines(tmp[,2], tmp[,4], lwd=3, lty=3,...)
    ##return(tmp)
}



### selecting 
idx <- 1:450
m.true <- as.numeric(colSums(df.noNoise))[idx]
dat <- df[,idx]

#### PCoA
dat.pcoa <- cmdscale(vegan::vegdist(t(dat)))
plot(dat.pcoa[,1], dat.pcoa[,2], col=c(rep(1, 400), rep(2, length(idx)-400)))


ncpu=10
max.iter <- 80
scaling <- median(m.true)
dev <- 3

res <- func.EM(dat, ncpu=ncpu, scaling=scaling, dev=dev, max.iter=max.iter, refine.start.iter = 30)

##res <- func.EM(dat, ncpu=ncpu, scaling=scaling, dev=4, max.iter=max.iter, refine.start.iter = 30)

trace.m <- res$trace.m
trace.p <- res$trace.p
err <- res$err
plot_d(trace.m, trace.p)


diag(b.truth) <- NA
est <- matrix(trace.p[-c(1:p),ncol(trace.p)], p,p)
diag(est) <- NA
par(mfrow = c(2,1))
auc.ad(est, b.truth)
pr.ad(est, b.truth)

pca <- prcomp( log(abs(err) + 1) )
biplot(pca, xlim=c(-0.5,0.5))

axis1 <- pca$x[,1]
axis2 <- pca$x[,2]

plot(axis1, axis2)
points(axis1[401:nrow(pca$x)], axis2[401:nrow(pca$x)], col='red', pch=19)

## center
points(median(axis1), median(axis2), pch=17, col='darkgreen', cex=3)
## axis1
abline(v=c(median(axis1) + 2*mad(axis1), median(axis1) - 2*mad(axis1)))
## axis2
abline(h=c(median(axis2) + 2*mad(axis2), median(axis2) - 2*mad(axis2)))

bad.idx <- rowSums(apply(pca$x[,1:2], 2, function(x) abs(x-median(x))/mad(x) > 2)) > 0

table(bad.idx)
table(bad.idx[401:450])

##############################################
## second run:
idx <- idx[!bad.idx]
## compare
idx <- idx[idx <= 400]
m.true <- as.numeric(colSums(df.noNoise))[idx]
scaling <- median(m.true)
dat <- df[,idx]

res <- func.EM(dat, ncpu=ncpu, scaling=scaling, dev=2, max.iter=max.iter)

trace.m <- res$trace.m
trace.p <- res$trace.p
err <- res$err
plot_d(trace.m, trace.p)

diag(b.truth) <- NA
est <- matrix(trace.p[-c(1:p),ncol(trace.p)], p,p)
diag(est) <- NA
par(mfrow = c(2,1))
auc.ad(est, b.truth)
pr.ad(est, b.truth)

pca <- prcomp(log(abs(err) + 1))

axis1 <- pca$x[,1]
axis2 <- pca$x[,2]

plot(axis1, axis2)
points(axis1[which(idx>400)], axis2[idx>400], col='red', pch=19)



## diagnostics:
par(mfrow = c(2,1))
boxplot(err)
rownames(err) <- NULL
colnames(err) <- NULL
pheatmap(t(log10(abs(err+1e-5))), cluster_rows = F, cluster_cols = F)

plot(rowSums(abs(err) > 0.1), type='l')
##matplot((err), cex=1)

####################### which points are wrong ##################
### noise free
dat.tss <- tss(df.noNoise)
i = 9
fil <- dat.tss[i,]!=0 #& !sample.filter[i,]
diag(b.truth) <- -1
tmp <- as.numeric(a.truth[i] + (b.truth[i,] %*% dat.tss[,fil])*m.true[fil])
plot(tmp)

##idx.true <- (tmp < -0.2)
idx.zero <- rep(FALSE, ncol(dat.tss))
idx.zero[dat.tss[i,]==0] <- TRUE
idx.zero[!idx.zero][tmp < -0.2] <- TRUE

df.fil <- df 
df.fil[9, idx.zero] <- 0
res <- func.EM(df.fil, ncpu=ncpu, scaling=scaling, dev=dev, max.iter=100)
plot_d(res$trace.m, res$trace.p)
#################################################################

tmp <- err[,9]
plot(abs(tmp))
points((1:length(tmp))[idx.zero], abs(tmp[idx.zero]), pch=19)

## redo:
df.fil <- df
index <- abs(err[,9])  <  1
points((1:length(tmp))[index], abs(tmp[index]), col='red')

df.fil[9, index] <- 0
res <- func.EM(df.fil, ncpu=ncpu, scaling=scaling, dev=dev, max.iter=max.iter)
plot_d(res$trace.m, res$trace.p)
##png("sp15.png", width=800, height=400, res=100)
##########################
## ad-hoc analysis with SPIEC-EAIS mb similar method

dat <- df
dat.tss <- tss(dat)

registerDoMC(10)

infer.ad <- function(Y, X, intercept=TRUE){
        fit <- cv.glmnet(X, Y, intercept=intercept)
        return(as.numeric(coef(fit))[-1])
}

res <- foreach(i=1:nrow(dat.tss), .combine=rbind) %dopar% {
    fil <- dat.tss[i,]!=0 
    X <- t(dat.tss[-i,][,fil])
    Y <- dat.tss[i, fil]
    theta <- rep(0, nrow(dat.tss))
    theta[i] <- -1 ## -beta_{ii}/beta_{ii}
    theta[-(i)] <- infer.ad(Y, X, intercept = TRUE)
    theta
}

library(SpiecEasi)

se.est.mb <- spiec.easi(t(dat.tss), method='mb', lambda.min.ratio=0.01, nlambda=20, icov.select.params=list(rep.num=50))
out.mb <- as.matrix(with(se.est.mb, (beta[[opt.index]] + t(beta[[opt.index]]))/2))
stability.mb <- as.matrix(se.est.mb$merge[[se.est.mb$opt.index]])

auc.ad <- function(est, lab, ...){
    est <- abs(est[!is.na(est)])
    lab <- lab[!is.na(lab)]
    lab <- (lab!=0 ) *1
    plot.roc(lab,est, print.auc = TRUE,...)
    pROC::auc(lab,est)
}


pr.ad <- function(esti, lab, add=FALSE, ...){
    ## lab=c(b.truth)
    ##esti = c(est.regress)
    esti <- abs(esti[!is.na(esti)])

    lab <- lab[!is.na(lab)]
    lab <- (lab!=0 ) *1
    ordered <- order(esti, decreasing = TRUE)
    ## wrong signs are not considered as tps
    m <- cbind(esti, lab)[ordered,]
    n_pos <- sum((m[,2]!=0) * 1)
    n_tp <- cumsum(m[,2])
    recall <- n_tp/n_pos
    precision <- n_tp/(1:NROW(m))
    f1 <- 2*(recall*precision)/(precision+recall)
    tmp <- cbind(thre=esti[ordered], recall, precision, f1)
    if(add==FALSE){
        plot(tmp[,2], tmp[,3], ylab='Precision',xlab='Recall', type='l',
             lwd=3, cex.axis=1.6, cex.lab=1.6,...)
    }else{
        lines(tmp[,2], tmp[,3], type='l',
             lwd=3, cex.axis=1.6, cex.lab=1.6,...)
    }
    lines(tmp[,2], tmp[,4], lwd=3, lty=3,...)
    ##return(tmp)
}


par(mfrow = c(1,3))
diag(b.truth) <- NA
est.regress <- as.matrix(res)
diag(est.regress) <- NA
plot(c(b.truth[est.regress<200]),c(est.regress[est.regress<200]) , main='Parameters (b)',
     xlab='Truth',ylab='Regression (MB-like)', cex=2, ##ylim=c(-5,8), xlim=c(-5,8),
     pch=19, col=rgb(0,0,0,0.3),cex.axis=1.5,cex.lab=1.5)
abline(0,1, lwd=2,lty=2)

est.se <- as.matrix(out.mb)
diag(est.se) <- NA
diag(stability.mb) <- NA
plot(c(b.truth),c(est.se) , main='Parameters (b)',
     xlab='Truth',ylab='SPIE-EASI (MB)', cex=2,
     pch=19, col=rgb(0,0,0,0.3),cex.axis=1.5,cex.lab=1.5)
abline(0,1, lwd=2,lty=2)

est.beem <- matrix(trace.p[-c(1:p),ncol(trace.p)], p,p)
diag(est.beem) <- NA
plot(c(b.truth),c(est.beem) , main='Final parameters (b)',
     xlab='Truth',ylab='BEEM', cex=2, ylim=c(-5,8), xlim=c(-5,8),
     pch=19, col=rgb(0,0,0,0.3),cex.axis=1.5,cex.lab=1.5)
abline(0,1, lwd=2,lty=2)


auc.ad(c(est.regress),c(b.truth), col='darkgreen')
auc.ad(c(est.beem), c(b.truth), col='black', add=TRUE, print.auc.y=0.8, print.auc.x=0.8)
auc.ad(c(est.se), c(b.truth), col='purple', add=TRUE, print.auc.x=0.8)

pr.ad(c(est.regress), c(b.truth), col='darkgreen', ylim=c(0,1))
pr.ad(c(est.beem), c(b.truth), col='black', add=TRUE)
pr.ad(c(stability.mb), c(b.truth), col='purple', add=TRUE)



plot(c(abs((est-b.truth)/b.truth)), ylim=c(0,1.2))

##dev.off()

par(cex=2, mfrow = c(1,2))

plot(a.truth, trace.p[1:p, 2], main='Final parameters (a)',
     xlab='Truth',ylab='BEEM', pch=19, col='red',cex=2
     )
abline(0,1, lwd=2,lty=2)

diag(b.truth) <- NA
est <- matrix(trace.p[-c(1:p), 2], p,p)
diag(est) <- NA

plot(c(b.truth),c(est) , main='Final parameters (b)',
     xlab='Truth',ylab='BEEM', pch=19, col='blue',cex=2)
abline(0,1, lwd=2,lty=2)

eval.func()

