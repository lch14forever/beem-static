##' @title resample.EM
##'
##' @param data data in the format of cbind(Y, X)
##' @param external.perturbation external perturbation presence matrix (each perturbation in one row, each sample in one column) (Default: NULL)
##' @param m biomass initialization
##' @param perc percentage of samples to take for each iteration
##' @param res.iter number of resample iteration
##' @param ... additional parameters for beemStatic::func.EM
##' @description (deprecated) Resampling for the inference process
##' @export
##' @author Chenhao Li, Gerald Tan, Niranjan Nagarajan
resample.EM <- function(data, external.perturbation = NULL, m, perc, res.iter, ...) {
  n <- ncol(data)
  p <- nrow(data)
  res_p <- data.frame(matrix(NA, p*(p+1), res.iter))
  res_m <- data.frame(matrix(NA, n, res.iter))

  for(i in 1:res.iter){
    message(paste0("#### Resample iteration: ", i, " #####"))
    indices <- sort(sample(1:n, n*perc))
    if (!is.null(external.perturbation)) {
      num.perturb <- nrow(external.perturbation)
      tmp <- func.EM(data[, indices], external.perturbation = external.perturbation[, indices], m.init=m[indices], ...)
    } else {
      tmp <- func.EM(data[, indices], external.perturbation = NULL, m.init=m[indices], ...)
    }
    res_p[,i] <- tmp$trace.p[, ncol(tmp$trace.p)]
    res_m[indices,i] <- tmp$trace.m[, ncol(tmp$trace.m)]
  }
  if (!is.null(external.perturbation)) {
    list(res.a = res_p[1:p,], res.b = res_p[(p+1):(nrow(res_p)-num.perturb*p),], res.c = res_p[(nrow(res_p)-num.perturb*p + 1):nrow(res_p),], res.m = res_m)
  } else {
    list(res.a = res_p[1:p,], res.b = res_p[-c(1:p),], res.m = res_m)
  }

  list(res.a = res_p[1:p,], res.b = res_p[-c(1:p),])
}
