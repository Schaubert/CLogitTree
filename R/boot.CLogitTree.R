#' Bootstrap confidence intervals for CLogitTree or CLogitForest
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
#' @param quant.type Algorithm type used for quantiles, as described in argument \code{type} from \code{\link{quantile}}.
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
#'                         perm_test = FALSE, depth_max=4)
#' illu.tree <- pruneBIC(illu.tree)
#'
#' set.seed(1860)
#' tree.boot <- bootci(illu.tree, B = 10)
#' tree.boot
#' }
#' @export
bootci.CLogitTree <- function(model, B = 100, ncores = 2, alpha.ci = 0.05, quant.type = 7, ...){

  Xboot = cbind(model$Z, model$exposure, model$strata, model$y)
  names(Xboot)[(ncol(Xboot)-2):ncol(Xboot)] <- c("Xexpo", "Xstrata","Xy")

  index <- levels(droplevels(as.factor(model$strata)))

  alpha <- model$alpha
  epsilon <- model$epsilon
  nperm <- model$nperm
  minnodesize <- model$minnodesize
  minbucket <- model$minbucket
  depth_max <- model$depth_max
  perm_test <- model$perm_test
  mtry <- model$mtry
  lambda <- model$lambda
  print.trace <- model$print.trace
  fit <- model$fit
  prunedBIC <- model$prunedBIC

  if(is.null(model$beta_hat)){
    mod <- CLogitTree(data = Xboot,
                       response = "Xy",
                       exposure="Xexpo",
                       s = "Xstrata",
                       alpha = alpha,
                       nperm = nperm,
                       minnodesize=minnodesize,
                       minbucket = minbucket,
                       depth_max = depth_max,
                       perm_test = perm_test,
                       mtry = mtry,
                       lambda = lambda,
                       print.trace = print.trace,
                       fit = TRUE,
                       ncores = ncores,
                       epsilon = epsilon)$beta_hat
    beta_hat <- mod$beta_hat

    if(prunedBIC){
      beta_hat <- pruneBIC(mod)$beta_hat
    }
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
  #                   fit=TRUE,
  #                   epsilon = epsilon)
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
                  minbucket = minbucket,
                  depth_max = depth_max,
                  perm_test = perm_test,
                  mtry = mtry,
                  lambda = lambda,
                  print.trace = print.trace,
                  fit = fit,
                  prunedBIC = prunedBIC,
                  epsilon = epsilon)
  }else{

    cl <- makeCluster(ncores, outfile = "")

    clusterExport(cl, varlist = c("index", "Xboot", "alpha", "nperm", "minnodesize",
                                  "minbucket", "depth_max", "perm_test", "mtry",
                                  "lambda", "print.trace", "fit", "prunedBIC", "epsilon"),
                  envir = sys.frame(sys.nframe()))

    beta.hat.vec <- parSapply(cl, seeds, one_boot_fun,
                     index = index,
                     Xboot = Xboot,
                     alpha = alpha,
                     nperm = nperm,
                     minnodesize = minnodesize,
                     minbucket = minbucket,
                     depth_max = depth_max,
                     perm_test = perm_test,
                     mtry = mtry,
                     lambda = lambda,
                     print.trace = print.trace,
                     fit = fit,
                     prunedBIC = prunedBIC,
                     epsilon = epsilon)
    stopCluster(cl)
  }


  return(list(beta_hat = beta_hat, ci_beta = quantile(beta.hat.vec, c(alpha.ci/2, 1-alpha.ci/2),
                                                      type = quant.type)))
}



