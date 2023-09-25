#' Bootstrap confidence intervals for CLogitForest
#'
#' Performs estimation of bootstrap confidence intervals for separate exposure variable effects from \code{\link{CLogitForest}} objects.
#'
#' @param model Original CLogitForest model
#' @param B Number of bootstrap iterations
#' @param ncores Number of parallel cores to use
#' @param alpha.ci Confidence level of confidence intervals is calculated as \code{1-alpha.ci}.
#' @param ... Further bootci arguments
#' @return
#' \item{beta_hat}{Estimate for separate exposure  effect}
#' \item{ci_beta}{Estimated bootstrap confidence interval for exposure effect}
#' @author Gunther Schauberger: \email{gunther.schauberger@@tum.de} \cr
#' Moritz Berger: \email{moritz.berger@@imbie.uni-bonn.de}
#' @seealso \code{\link{CLogitForest}}
#' @examples
#' \donttest{
#' data(illu.small)
#'
#' set.seed(1860)
#' illu.rf <- CLogitForest(illu.small, response = "y", exposure = "x", s = "strata",
#'                         ntree = 4, depth_max=2, tune.mtry = FALSE)
#'
#' set.seed(1860)
#' rf.boot <- bootci(illu.rf, B = 10)
#' rf.boot
#' }
#' @export
bootci.CLogitForest <- function(model, B = 100, ncores = 1, alpha.ci = 0.05, ...){


  if(is.null(model$exposure)){
    stop("Confidence intervals only make sense for model where an explicit exposure effect is specified!")
  }

# browser()
  ntree <- model$ntree
  epsilon <- model$epsilon
  mtry <- model$mtry
  subsample.632 <- model$subsample.632
  minnodesize <- model$minnodesize
  minbucket <- model$minbucket
  depth_max <- model$depth_max
  lambda <- model$lambda
  print.trace <- model$print.trace
  fit <- model$fit
  response <- model$response
  s <- model$s
  exposure <- model$exposure
  data <- model$data
  BIC <- model$BIC
  # offset <- model$offset
  linear.offset <- model$linear.offset


  index <- levels(droplevels(as.factor(data[,names(data)==s])))


  if(is.null(model$beta_hat)){
    mod <- CLogitForest(data = data,
                       response = response,
                       exposure= exposure,
                       s = s,
                       ntree = ntree,
                       mtry = mtry,
                       subsample.632 = subsample.632,
                       minnodesize=minnodesize,
                       minbucket = minbucket,
                       depth_max = depth_max,
                       lambda = lambda,
                       print.trace = print.trace,
                       fit = TRUE,
                       ncores = ncores,
                       BIC = BIC,
                       #offset = offset,
                       linear.offset = linear.offset,
                       epsilon = epsilon)$beta_hat
  }else{
    beta_hat <- model$beta_hat
  }


  seeds <- abs(round(rnorm(B) * 1e8))


  if(ncores ==1){
    beta.hat.vec <- sapply(seeds, one_boot_forest,
                  index = index,
                  Xboot = data,
                  minnodesize = minnodesize,
                  minbucket = minbucket,
                  depth_max = depth_max,
                  ntree = ntree,
                  mtry = mtry,
                  subsample.632 = subsample.632,
                  lambda = lambda,
                  print.trace = print.trace,
                  fit = fit,
                  response = response,
                  exposure= exposure,
                  s = s,
                  BIC = BIC,
                  # offset = offset,
                  linear.offset = linear.offset,
                  epsilon = epsilon)
  }else{

    # cl <- makeCluster(ncores, outfile = "CLogitForest_bootlog.txt")
    cl <- makeCluster(ncores, outfile = "")

    clusterExport(cl, varlist = c("index", "data", "minnodesize", "ntree",
                                  "minbucket", "depth_max", "mtry",
                                  "lambda", "print.trace", "fit", "exposure",
                                  "response", "s", "BIC", "linear.offset","epsilon"),
                  envir = sys.frame(sys.nframe()))
    clusterEvalQ(cl, library(CLogitTree))

    beta.hat.vec <- parSapply(cl, seeds, one_boot_forest,
                              index = index,
                              Xboot = data,
                              minnodesize = minnodesize,
                              minbucket = minbucket,
                              depth_max = depth_max,
                              ntree = ntree,
                              mtry = mtry,
                              subsample.632 = subsample.632,
                              lambda = lambda,
                              print.trace = print.trace,
                              fit = fit,
                              response = response,
                              exposure= exposure,
                              s = s,
                              BIC = BIC,
                              # offset = offset,
                              linear.offset = linear.offset,
                              epsilon = epsilon)
    stopCluster(cl)
  }


  return(list(beta_hat = beta_hat, ci_beta = quantile(beta.hat.vec, c(alpha.ci/2, 1-alpha.ci/2)),
              beta_hat_vector = beta.hat.vec))
}



