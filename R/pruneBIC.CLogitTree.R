#' Pruning for CLogitTree
#'
#' @param tree Original CLogitTree model
#' @param ... Further arguments
pruneBIC <- function (tree, ...) {
  UseMethod("pruneBIC", tree)
}
#'
#' Pruning for CLogitTree
#'
#' Prunes trees built by \code{\link{CLogitTree}}. Pruning is steered by the BIC of the sequence of fitted models during tree building.
#'
#' @param tree Original CLogitTree model
#' @param ... Further prune arguments
#' @return \code{\link{CLogitTree}} object
#' @author Gunther Schauberger: \email{gunther.schauberger@@tum.de} \cr
#' Moritz Berger: \email{moritz.berger@@imbie.uni-bonn.de}
#' @seealso  \code{\link{CLogitTree}}, \code{\link{plot.CLogitTree}}, \code{\link{bootci.CLogitTree}}
#' @examples
#' data(illu.small)
#'
#' set.seed(1860)
#' illu.tree <- CLogitTree(illu.small, response = "y", exposure = "x", s = "strata",
#'                         perm_test = FALSE, print.trace = FALSE, depth_max = 2)
#'
#' plot(illu.tree)
#' plot(pruneBIC(illu.tree))
#'
#' @export
pruneBIC.CLogitTree <- function(tree, ...){

  iter <- which.min(tree$BIC)

  model      <- tree$model[[iter]]
  params     <- tree$param[[iter]]
  params_fit <- tree$param_fit[[iter]]
  design     <- tree$design[[iter]]
  pvalues    <- tree$pvalue[1:iter]
  devs       <- tree$dev[1:iter]
  crits      <- tree$crit[1:iter]
  BICs       <- tree$BIC[1:iter]
  y_tab      <- tree$y_tab[1:iter]

  if(iter==1){
    if(length(tree$exposure)>0){
      beta_hat <- unlist(model$penalized)
    } else{
      beta_hat <- NULL
    }
    gamma_hat <- NULL
    splits   <- NULL
  }

  if(iter>1){
    if(length(tree$exposure)>0){
      gamma_hat <- c(unlist(model$penalized),0)
      beta_hat  <- unlist(model$unpenalized)
      names(gamma_hat) <- params_fit[-1]
      names(beta_hat)  <- params_fit[1]
    } else{
      gamma_hat <- c(unlist(model$penalized),0)
      beta_hat  <- NULL
      names(gamma_hat) <- params_fit
    }
    gamma_hat <- gamma_hat[params]
    splits <- tree$splits[1:(iter-1),]
  }

  if(iter>1){
    gamma_hat_sym <- gamma_hat - mean(gamma_hat)
  } else{
    gamma_hat_sym <- NULL
  }

  to_return <-  list("beta_hat"=beta_hat,
                     "gamma_hat"=gamma_hat,
                     "gamma_hat_sym" = gamma_hat_sym,
                     "splits"=splits,
                     "Z"=tree$Z,
                     "y"=tree$y,
                     "y_tab"=y_tab,
                     "strata" = tree$strata,
                     "exposure" = tree$exposure,
                     "model"=model,
                     "design"=design,
                     "param"=params,
                     "param_fit"=params_fit,
                     "pvalue"=pvalues,
                     "dev"=devs,
                     "crit"=crits,
                     "BIC"=BICs,
                     "minnodesize" = tree$minnodesize,
                     "minbucket" = tree$minbucket,
                     "depth_max" = tree$depth_max,
                     "lambda" = tree$lambda,
                     "trace" = tree$trace,
                     "nperm" = tree$nperm,
                     "alpha" = tree$alpha,
                     "perm_test" = tree$perm_test,
                     "mtry" = tree$mtry,
                     "call"=tree$call,
                     "print.trace" = tree$print.trace,
                     "prunedBIC" = TRUE)

  class(to_return) <- "CLogitTree"
  return(to_return)
}
