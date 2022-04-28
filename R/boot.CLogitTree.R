#' Bootstrap confidence intervals for CLogitTree
#'
#' @param model Object to calculate bootstrap confidence intervals for
#' @param ... Further arguments
bootci <- function (model, ...) {
  UseMethod("bootci", model)
}
#' Bootstrap confidence intervals for CLogitTree
#'
#' Performs estimation of bootstrap confidence intervals for separate exposure variable effects from \code{\link{CLogitTree}} objects.
#'
#' @param model Original CLogitTree model
#' @param B Number of bootstrap iterations
#' @param ncores Number of parallel cores to use
#' @param alpha.ci Confidence level of confidence intervals is calculated as \code{1-alpha.ci}.
#' @param ... Further bootci arguments
#' @return
#' \item{beta_hat}{Estimate for separate exposure  effect}
#' \item{ci_beta}{Estimated bootstrap confidence interval for exposure effect}
#' @author Gunther Schauberger: \email{gunther.schauberger@@tum.de} \cr
#' Moritz Berger: \email{moritz.berger@@imbie.uni-bonn.de}
#' @seealso \code{\link{CLogitTree}}, \code{\link{plot.CLogitTree}}, \code{\link{prune.CLogitTree}}
#' @examples
#' \donttest{
#' data(illu.small)
#'
#' set.seed(1860)
#' illu.tree <- CLogitTree(illu.small, response = "y", exposure = "x", s = "strata",
#'                         alpha = 0.05, nperm = 20, print.trace = FALSE)
#'
#' plot(illu.tree)
#'
#' set.seed(1860)
#' illu.boot <- bootci(illu.tree, B = 10, alpha = 0.2)
#' }
#' @export
bootci.CLogitTree <- function(model, B = 100, ncores = 2, alpha.ci = 0.05, ...){

  Xboot = cbind(model$Z, model$exposure, model$strata, model$y)
  names(Xboot)[(ncol(Xboot)-2):ncol(Xboot)] <- c("Xexpo", "Xstrata","Xy")

  index <- levels(droplevels(as.factor(model$strata)))

  alpha <- model$alpha
  nperm <- model$nperm
  minnodesize <- model$minnodesize
  perm_test <- model$perm_test
  mtry <- model$mtry
  lambda <- model$lambda
  print.trace <- model$print.trace
  fit <- model$fit

  if(is.null(model$beta_hat)){
    beta_hat <- CLogitTree(data = Xboot,
                       response = "Xy",
                       exposure="Xexpo",
                       s = "Xstrata",
                       alpha = alpha,
                       nperm = nperm,
                       minnodesize=minnodesize,
                       perm_test=perm_test,
                       mtry=mtry,
                       lambda=lambda,
                       print.trace=print.trace,
                       fit=TRUE,
                       ncores = ncores)$beta_hat
  }else{
    beta_hat <- model$beta_hat
  }



  # beta.hat.vec <- rep(NA,B)
  #
  # for(j in 1:B){
  #   cat("Starting bootstrap iteration ",j,"out of ",B,"\n")
  # index2 <- sample(index, size = length(index), replace = TRUE)
  #
  #
  #
  # X2 <- c()
  # for(i in 1:length(index2)){
  #   Xnew <- X[X$Xstrata == index2[i],]
  #   Xnew$Xstrata <- rep(i, nrow(Xnew))
  #   X2 <- rbind(X2, Xnew)
  # }
  #
  # ret <- CLogitTree(data = X2,
  #                   response = "Xy",
  #                   exposure="Xexpo",
  #                   s = "Xstrata",
  #                   alpha = model$alpha,
  #                   nperm = model$nperm,
  #                   minnodesize=model$minnodesize,
  #                   perm_test=model$perm_test,
  #                   mtry=model$mtry,
  #                   lambda=model$lambda,
  #                   print.trace=model$print.trace,
  #                   fit=TRUE)
  #
  # beta.hat.vec[j] <- ret$beta_hat
  # }



  seeds <- abs(round(rnorm(B) * 1e8))

  if(ncores ==1){
    beta.hat.vec <- sapply(seeds, one_boot_fun,
                  index = index,
                  Xboot = Xboot,
                  alpha = alpha,
                  nperm = nperm,
                  minnodesize = minnodesize,
                  perm_test = perm_test,
                  mtry = mtry,
                  lambda = lambda,
                  print.trace = print.trace,
                  fit = fit)
  }else{

    cl <- makeCluster(ncores, outfile = "")

    clusterExport(cl, varlist = c("index", "Xboot", "alpha", "nperm", "minnodesize",
                                  "perm_test", "mtry", "lambda", "print.trace", "fit"),
                  envir = sys.frame(sys.nframe()))

    beta.hat.vec <- parSapply(cl, seeds, one_boot_fun,
                     index = index,
                     Xboot = Xboot,
                     alpha = alpha,
                     nperm = nperm,
                     minnodesize = minnodesize,
                     perm_test = perm_test,
                     mtry = mtry,
                     lambda = lambda,
                     print.trace = print.trace,
                     fit = fit)
    stopCluster(cl)
  }



  return(list(beta_hat = beta_hat, ci_beta = quantile(beta.hat.vec, c(alpha.ci/2, 1-alpha.ci/2))))
}






# boot.CLogitTree(data, model, R){
#
#   boot(data, boot.CLfun, R = R, model = model)
# }
#
# set.seed(1860)
# b1 <- boot.CLogitTree2(sub.mat2, tree.mod2)
#
# set.seed(1860)
# boot::boot(sub.mat2, boot.CLfun, R = 10, model = tree.mod2)
