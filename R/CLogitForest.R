#' Function to perform CLogitForest
#'
#' Performs CLogitForest, a random forest method for the analysis of matched case-control studies. The estimation is embedded into the framework of conditional logistic regression.
#' Splits are incorporated into the model as separate dummy variables specifying the terminal nodes. Additional methods for pruning, plotting and bootstrap are available.
#'
#' @param data Data frame (containing all required variable)
#' @param response Name of response variable (character)
#' @param exposure Name of (optional) separate exposure variable (character)
#' @param s Name of strata variable (character)
#' @param ntree Number of trees to fit for random forest
#' @param mtry Number of variables that a randomly selected for potential splits. Only effective if \code{tune.mtry = FALSE}
#' @param subsample.632 If \code{TRUE}, only a subsample of 63.2% of all strata is
#' randomly selected for a single tree
#' @param minnodesize Minimal node size in order to be eligible for further splitting
#' @param minbucket Minimum number of observations in any terminal node
#' @param depth_max Maximum depth of the tree, with the root node counted as depth 0. If \code{NULL} (default), the size of the trees is not restricted.
#' @param lambda Tuning parameter for optional L2 penalty
#' @param print.trace Shall trace of permutation tests be printed?
#' @param ncores Number of parallel nodes to use.
#' @param fit Shall the internally fitted models be returned (required for the use of \code{\link[CLogitTree]{prune}})?
#' @param BIC Shall each tree be pruned using the BIC criterion?
#' @param linear.offset Shall a linear offset (from regular CLR) be used?
#' @param tune.mtry Shall \code{mtry} be tuned prior to model fitting? If true, \code{tune.ntree} trees are fitted with all possible
#' values of \code{mtry}
#' @param tune.ntree Number of trees to be used for \code{mtry} tuning. Only effective if \code{tune.mtry = TRUE}.
#' @param tune.range Range of possible values for mtry to be considered mtry-tuning. Only effective if \code{tune.mtry = TRUE}.
#' @param epsilon Convergence tolerance. Iteration continues until the relative change
#' in the conditional log likelihood is less than eps. Must be positive.
#'
#' @return
#' \item{beta_hat}{Estimate for separate exposure  effect (arithmetic mean across exposure effects of all trees)}
#' \item{beta_median_hat}{Median across exposure effects of all trees}
#' \item{beta_list}{Vector containing all exposure estimates}
#' \item{tree_list}{List containing all CLogitTree objects which constitute the tree}
#' \item{mtry}{mtry argument from call}
#' \item{call}{function call}
#' \item{ncores}{Number of parallel nodes}
#' \item{minnodesize}{Minimal node size in order to be eligible for further splitting}
#' \item{minbucket}{Minimum number of observations in any terminal node}
#' \item{depth_max}{Maximum depth of the tree}
#' \item{lambda}{L2 tuning parameter}
#' \item{print.trace}{trace argument from call}
#' \item{subsample.632}{If \code{TRUE}, only a subsample of 63.2% of all strata was
#' randomly selected for a single tree}
#' \item{data}{Original data set}
#' \item{response}{Name of response variable}
#' \item{exposure}{Name of exposure variable}
#' \item{s}{Name of strata variable}
#' \item{ntree}{Number of trees to fit for random forest}
#' \item{BIC}{BIC argument from call}
#' \item{lin.off.model}{CLR model fitted for linear offset (if \code{linear.offset = TRUE})}
#' \item{offset.model}{If \code{linear.offset = TRUE}, this is the finally fitted model with the random forest as its offset.}
#' \item{epsilon}{Convergence tolerance}
#'
#' @author Gunther Schauberger: \email{gunther.schauberger@@tum.de} \cr
#' Moritz Berger: \email{Moritz.Berger@@imbie.uni-bonn.de}
#' @seealso \code{\link{bootci.CLogitForest}}, \code{\link{predict.CLogitForest}}, \code{\link{varimp.CLogitForest}}
#' @examples
#' data(illu.small)
#'
#' set.seed(1860)
#' illu.rf <- CLogitForest(illu.small, response = "y", exposure = "x", s = "strata",
#'                         ntree = 4, depth_max=2, tune.mtry = FALSE)
#'
#' illu.rf
#' @export
CLogitForest <- function(data,
                       response,
                       exposure = NULL,
                       s,
                       ntree = 100,
                       mtry = ifelse(is.null(exposure),max(floor((ncol(data)-2)/3), 1),
                                     max(floor((ncol(data)-3)/3), 1)),
                       subsample.632 = FALSE,
                       minnodesize = 15,
                       minbucket = 5,
                       depth_max = NULL,
                       lambda = 1e-20,
                       print.trace = FALSE,
                       fit = TRUE,
                       ncores = 2,
                       BIC = FALSE,
                       linear.offset = FALSE,
                       tune.mtry = TRUE,
                       tune.ntree = 20,
                       tune.range = 2:floor((ncol(data)-2)*0.7),
                       epsilon = 1e-5){

  ## preliminary definitions
  n    <- nrow(data)
  nvar <- ncol(data[,!names(data)%in%c(exposure,s)])

  seeds <- abs(round(rnorm(ntree) * 1e8))


if(linear.offset){

  names.Z    <- names(data)[!(names(data)%in%c(response, exposure,s))]

  if(is.null(exposure)){
    form.linear <- as.formula(paste0(response," ~ ", paste0(names.Z,collapse="+")))
  }else{
    form.linear <- as.formula(paste0(response," ~ ", exposure, "+", paste0(names.Z,collapse="+")))
  }

  Z <- model.matrix(form.linear, data = data)[,-1,drop = FALSE]

  if(!is.null(exposure)){
    invisible(capture.output(lin.off.model <- penalized.clr(response=data[,response], stratum=data[,s],
                                   penalized=Z[,-1,drop = FALSE], unpenalized=Z[,1],
                                   lambda=lambda, alpha=10e-20, epsilon = epsilon)))

    linear.off <- Z[,-1,drop = FALSE]%*%lin.off.model$penalized
  } else{
    invisible(capture.output(lin.off.model <- penalized.clr(response=data[,response], stratum=data[,s],
                                   penalized=Z, unpenalized=NULL,
                                   lambda=lambda, alpha=10e-20, epsilon = epsilon)))

    linear.off <- Z%*%lin.off.model$penalized
  }


}else{
  linear.off <- rep(0, nrow(data))
  lin.off.model <- NA
}

  if(tune.mtry){
    mtry <- tune_mtry(seeds, subsample.632, data,
                      s, exposure, response, minnodesize, minbucket,
                      depth_max, lambda, print.trace, fit, BIC, ncores, tune.ntree, tune.range,
                      linear.off, epsilon)
  }
  # browser()

  ## one core or several cores
  if(ncores == 1){
    tree_list <- lapply(seeds, clrf.fit,
                         subsample.632 = subsample.632,
                         data = data,
                         s = s,
                         exposure = exposure,
                         response = response,
                         minnodesize = minnodesize,
                         minbucket = minbucket,
                         depth_max = depth_max,
                         mtry = mtry,
                         lambda = lambda,
                         print.trace = print.trace,
                         fit = fit,
                         BIC = BIC,
                         linear.off = linear.off,
                        epsilon = epsilon)
  }else{
    # cl <- makeCluster(ncores, outfile = "CLogitForest_log.txt")
    cl <- makeCluster(ncores, outfile = "")
    clusterExport(cl, varlist = c("data",
                         "response", "exposure", "s", "mtry", "subsample.632",
                         "minnodesize", "minbucket", "depth_max", "lambda",
                         "print.trace", "fit", "BIC", "linear.off", "epsilon"),
                  envir = sys.frame(sys.nframe()))

    tree_list <- parLapply(cl, seeds, clrf.fit,
                            subsample.632 = subsample.632,
                            data = data,
                            s = s,
                            exposure = exposure,
                            response = response,
                            minnodesize = minnodesize,
                            minbucket = minbucket,
                            depth_max = depth_max,
                            mtry = mtry,
                            lambda = lambda,
                            print.trace = print.trace,
                            fit = fit,
                            BIC = BIC,
                            linear.off = linear.off,
                           epsilon = epsilon)
    stopCluster(cl)
  }

  if(is.null(exposure)){
    exp.effect <- NA
    beta_hat <- NA
    beta_median_hat <- NA
  }else{
    exp.effect <- unlist(lapply(tree_list, function(x){x$tree.b$beta_hat}))


    beta_hat <- mean(exp.effect)
    beta_median_hat <- median(exp.effect)
  }


  to_return <- list("beta_hat"=beta_hat,
                    "beta_median_hat" = beta_median_hat,
                     "beta_list" = exp.effect,
                    "tree_list" = tree_list,
                     "call"=match.call(),
                     "mtry" = mtry,
                     "ncores" = ncores,
                     "minnodesize" = minnodesize,
                     "minbucket" = minbucket,
                     "depth_max" = depth_max,
                     "lambda" = lambda,
                     "print.trace" = print.trace,
                     "subsample.632" = subsample.632,
                     "data" = data,
                     "response" = response,
                     "exposure" = exposure,
                     "s" = s,
                     "ntree" = ntree,
                     "BIC" = BIC,
                      "lin.off.model" = lin.off.model,
                     "linear.offset" = linear.offset,
                      "epsilon" = epsilon)

  class(to_return) <- "CLogitForest"
  return(to_return)

}

clrf.fit <- function(seed, subsample.632, data,
                     s, exposure, response, minnodesize, minbucket,
                     depth_max, mtry, lambda, print.trace, fit, BIC, linear.off, epsilon){

  strata <- unique(data[, names(data) == s])

  n.strata <- length(strata)

  set.seed(seed)
  if(mtry == "random"){
    n.var <- ncol(data)-2
    if(!is.null(exposure)){
      n.var <- n.var-1
    }
    mtry <- sample(1:n.var,1)
  }

  # sample bootstrap sample or subsample
    if(subsample.632){
      sample.strata <- sample(strata, ceiling(0.632*n.strata), replace = FALSE)
      subset.index <- which(data[, names(data) == s] %in% sample.strata)
      data.sample <- data[subset.index,]
      linear.off <- linear.off[subset.index]
    }else{
      sample.strata <- sample(strata, n.strata, replace = TRUE)
      data.sample <- data[which(data[, names(data) == s] %in% sample.strata[1]),]
      data.sample[,s] <- 1
      off.help <- linear.off[which(data[, names(data) == s] %in% sample.strata[1])]
      for(iii in seq_along(sample.strata[-1])){
        data.add <- data[which(data[, names(data) == s] %in% sample.strata[iii+1]),]
        data.add[,s] <- iii+1
        data.sample <- rbind(data.sample, data.add)
        off.help <- c(off.help, linear.off[which(data[, names(data) == s] %in% sample.strata[iii+1])])
      }
      linear.off <- off.help
    }

  # index of indicating the selected rows of the data for the training data


    tree.b <- CLogitTree(data = data.sample,
                         response = response,
                         exposure = exposure,
                         s = s,
                         alpha = 0,
                         nperm = 1000,
                         minnodesize = minnodesize,
                         minbucket = minbucket,
                         depth_max = depth_max,
                         perm_test = FALSE,
                         mtry = mtry,
                         lambda = lambda,
                         print.trace = print.trace,
                         fit = fit,
                         offset = linear.off,
                         epsilon = epsilon)

    if(BIC){
      tree.b <- pruneBIC(tree.b)
    }

    return(list(tree.b = tree.b, sample.strata = sample.strata))
}


tune_mtry <- function(seeds, subsample.632, data,
                      s, exposure, response, minnodesize, minbucket,
                      depth_max, lambda, print.trace, fit, BIC, ncores, tune.ntree, tune.range,
                      linear.off, epsilon){

if(tune.ntree>length(seeds)){
  tune.ntree <- length(seeds)
}

n.var <- ncol(data)-2
if(!is.null(exposure)){
  n.var <- n.var-1
}

seeds <- seeds[1:tune.ntree]

  ## one core or several cores
  if(ncores == 1){
    ll_list <- sapply(seeds, tune_mtry_single,
                        subsample.632 = subsample.632,
                        data = data,
                        s = s,
                        exposure = exposure,
                        response = response,
                        minnodesize = minnodesize,
                        minbucket = minbucket,
                        depth_max = depth_max,
                        lambda = lambda,
                        print.trace = print.trace,
                        fit = fit,
                        BIC = BIC,
                        n.var = n.var,
                      tune.range = tune.range,
                        linear.off = linear.off,
                      epsilon = epsilon)
  }else{
    # cl <- makeCluster(ncores, outfile = "tunemtry_log.txt")
    cl <- makeCluster(ncores, outfile = "")

    clusterExport(cl, varlist = c("data",
                                  "response", "exposure", "s", "subsample.632",
                                  "minnodesize", "minbucket", "depth_max", "lambda",
                                  "print.trace", "fit","BIC","n.var", "linear.off",
                                  "epsilon","tune.range"),
                  envir = sys.frame(sys.nframe()))

    ll_list <- parSapply(cl, seeds, tune_mtry_single,
                           subsample.632 = subsample.632,
                           data = data,
                           s = s,
                           exposure = exposure,
                           response = response,
                           minnodesize = minnodesize,
                           minbucket = minbucket,
                           depth_max = depth_max,
                           lambda = lambda,
                           print.trace = print.trace,
                           fit = fit,
                           BIC = BIC,
                           n.var = n.var,
                         tune.range = tune.range,
                           linear.off = linear.off,
                         epsilon = epsilon)
    stopCluster(cl)
  }

ll <- rowSums(ll_list)

mtry.opt <- tune.range[which.max(ll)]
cat("The optimal mtry is ",mtry.opt,"\n")
mtry.opt
}


tune_mtry_single <- function(seed, subsample.632, data,
                      s, exposure, response, minnodesize, minbucket,
                      depth_max, lambda, print.trace, fit, BIC, n.var,
                      tune.range, linear.off, epsilon){



  ll <- rep(0,length(tune.range))

  strata <- unique(data[, names(data) == s])
  n.strata <- length(strata)

  sample.train <- sample(strata, ceiling(0.75*n.strata), replace = FALSE)
  subset.train <- which(data[, names(data) == s] %in% sample.train)
  data.train <- data[subset.train,]
  linear.off.train <- linear.off[subset.train]
  subset.test <- which(!(data[, names(data) == s] %in% sample.train))
  data.test <- data[subset.test,]


  for(tt in seq_along(tune.range)){
    fit.tt <- clrf.fit(seed, subsample.632, data.train,
                       s, exposure, response, minnodesize, minbucket,
                       depth_max, tune.range[tt], lambda, print.trace, fit, BIC, linear.off.train, epsilon)

    ll[tt] <- mean(predict(fit.tt$tree.b, type = "loglik",newdata = data.test,
                           response = response, exposure = exposure, s = s))
  }
  ll
}
