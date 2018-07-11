#!/mnt/software/bin/Rscript-3.1.2
args <- commandArgs(trailingOnly=T)
suppressMessages(library(pROC))
##library(pracma)
suppressMessages(library(reshape2))

rmse <- function(x, y){
    ## NOTE: multivariate RMSE only make sense when all variables are standardized (here everything is abundance)    
    sqrt(mean((scale(x) - scale(y))^2))
}

pr.mat <- function(response, predictor, descreasing = T){
    ordered <- order(predictor, decreasing = descreasing)
    ## wrong signs are not considered as tps
    m <- cbind(predictor, response)[ordered,]
    n_pos <- sum((m[,2]!=0) * 1)
    n_tp <- cumsum(m[,2])
    recall <- n_tp/n_pos
    precision <- n_tp/(1:NROW(m))
    f1 <- 2*(recall*precision)/(precision+recall)
    cbind(thre=predictor[ordered], recall, precision, f1)
}

rel_err <- function(x,y, pseudo=T){
    ##  x: truth
    ##  y: est
    ##unit <- min(abs(c(x[x!=0], y[y!=0])))
    ## pseudo <- 1e-9 ##10^floor(log10(unit))
    ## if(pseudo){
    ##     mean(abs((x-y)/(x+pseudo)))
    ## }else{
    ##     mean(abs((x-y)/(x)))}
    median(abs((x[x!=0]-y[x!=0])/(x[x!=0])))
}

eval.func <- function(parms.est, parms.truth, p, thre=1e-7){
    ## growth
    growth <- data.frame(Truth=parms.truth[1:p,3],
                         Est=parms.est[1:p,4])
    g.pearson <- cor(growth$Truth, growth$Est)
    g.spearman <- cor(growth$Truth, growth$Est, method='spearman')
    g.rmse <- rmse(growth$Truth, growth$Est)
    g.rel_err <- rel_err(growth$Truth, growth$Est)
    
    ## interactions
    interaction <- data.frame(Truth=parms.truth[(p+1):nrow(parms.truth),3],
                              Est=parms.est[(p+1):nrow(parms.truth),]$value)
    interaction.sig <- data.frame(Truth=parms.truth[(p+1):nrow(parms.truth),3],
                              Est=parms.est[(p+1):nrow(parms.truth),]$significance)
    ## self interations
    self.idx <- (parms.truth$Var2 == parms.truth$Var1)[-(1:p)]
    self.int <- interaction[self.idx,]
    s.pearson <- cor(self.int$Truth, self.int$Est)
    s.spearman <- cor(self.int$Truth, self.int$Est, method = 'spearman')
    s.rmse <- rmse(self.int$Truth, self.int$Est)
    s.rel_err <- rel_err(self.int$Truth, self.int$Est)
    a.rel_err <- rel_err(interaction$Truth, interaction$Est)
    ## interaction
    nonself.int <- interaction[!self.idx, ]
    nonself.int.sig <- interaction.sig[!self.idx, ]
    i.pearson <- cor(nonself.int$Truth, nonself.int$Est)
    i.spearman <- cor(nonself.int$Truth, nonself.int$Est, method='spearman')
    i.rmse <- rmse(nonself.int$Truth, nonself.int$Est)
    i.rel_err <- rel_err(nonself.int$Truth, nonself.int$Est)
    if(is.na(sum(interaction.sig))){
        i.auc <- as.numeric(auc(nonself.int$Truth!=0,
                            abs(nonself.int$Est)))
        i.pr.mat <- pr.mat(nonself.int$Truth!=0,
                       abs(nonself.int$Est))
    }else{
        i.auc <- as.numeric(auc(nonself.int$Truth!=0,
                                nonself.int.sig$Est))
        i.pr.mat <- pr.mat(nonself.int$Truth!=0,
                           nonself.int.sig$Est)
    }
    ##i.aupr <- trapz(i.pr.mat[,2], i.pr.mat[,3])
    tmp <- i.pr.mat[i.pr.mat[,1]>thre,]
    if(is.null(nrow(tmp))) {
        i.pr <- tmp
    }else{
        i.pr <- tail(tmp, 1)
    }
    if(NROW(i.pr)==0){
        i.precision <- i.recall <- i.f1 <- 0
    }else{
        i.precision <- i.pr[,3]
        i.recall <- i.pr[,2]
        i.f1 <- i.pr[,4]
    } 
    return(
        cbind(g.pearson,
              g.spearman,
              g.rel_err,
              s.pearson,
              s.spearman,
              s.rel_err,
              i.pearson,
              i.spearman,
              i.rel_err,
              i.auc,
              i.precision,
              i.recall,
              i.f1,
              a.rel_err
              )
        )
}

