##' @title auc.b
##' 
##' @param b.est
##' @param b.true
##' @description plot ROC curve with AUC for the interaction graph (interaction signs are ignored)
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
