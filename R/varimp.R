#' Variable importance for CLogitForest
#'
#' @param model Object to calculate variable importance for
#' @param ... Further arguments
varimp <- function (model, ...) {
  UseMethod("varimp", model)
}
#' Variable importance for CLogitForest
#'
#' Estimates variable importance for objects created by \code{\link{CLogitForest}}.
#'
#' @param model Original CLogitForest model
#' @param n.random Number of randomization iterations used per variable
#' @param oob Shall calculation of variable importance values only be based on out-of-bag observations
#' for each tree
#' @param ... Further varimp arguments
#' @return Estimated values of variable importance, separately for all variables
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
#' illu.rf
#'
#' vi <- varimp(illu.rf)
#' plot(vi)
#' }
#' @export
varimp.CLogitForest <- function(model, n.random = 1, oob = TRUE, ...){
  # browser()
## ncores option was taken out
  ncores <- 1


  # indicator if linear offset is used within RF
  linear.offset <- model$linear.offset
  off.model <- NA
  if(linear.offset){
    off.model <- model$lin.off.model
  }

  var.names <- names(model$data)[!(names(model$data)%in%c(model$response, model$exposure, model$s))]
  var.imp <- rep(0, length(var.names))
  names(var.imp) <- var.names

  var.list <- as.list(1:length(var.names))

if(oob){
# browser()
  single.trees <- lapply(model$tree_list, varimp.help.oob, s = model$s, data = model$data,
                       response = model$response, exposure = model$exposure,
                       linear.offset = linear.offset, off.model = off.model)

 if(ncores == 1){
   var.imp.list <- lapply(var.list, singlevar.oob,
                          model = model, var.names = var.names, single.trees = single.trees,
                          n.random = n.random,  linear.offset = linear.offset,
                          off.model = off.model)
 }else{
   # cl <- makeCluster(ncores, outfile = "CLogitForest_varimp_log.txt")
   cl <- makeCluster(ncores, outfile = "")

   clusterExport(cl, varlist = c("single.trees","var.list", "model", "var.names","n.random",
                                 "linear.offset", "off.model"),
                 envir = sys.frame(sys.nframe()))

   var.imp.list <- parLapply(cl, var.list, singlevar.oob,
                             model = model, var.names = var.names, single.trees = single.trees,
                             n.random = n.random,  linear.offset = linear.offset,
                             off.model = off.model)
   stopCluster(cl)

 }

 var.imp <- var.imp + unlist(var.imp.list)

}else{
  # browser()
  strata <- model$data[, names(model$data) == model$s]
  preds.orig <- predict(model, newdata = model$data, response = model$response,
                        exposure = model$exposure, s = model$s)
  preds.orig <- split(preds.orig, strata)
  y <- split(model$data[, names(model$data) == model$response], strata)


    if(ncores == 1){
      var.imp.list <- lapply(var.list, singlevar.ib,
                             model = model, var.names = var.names, preds.orig = preds.orig,
                             strata = strata, y=y, n.random = n.random)
    }else{
      cl <- makeCluster(ncores, outfile = "CLogitForest_varimp_log.txt")

      clusterExport(cl, varlist = c("preds.orig","var.list", "model", "var.names", "strata","y",
                                    "n.random"),
                    envir = sys.frame(sys.nframe()))

      var.imp.list <- parLapply(cl, var.list, singlevar.ib,
                                model = model, var.names = var.names, preds.orig = preds.orig,
                                strata = strata, y=y, n.random = n.random)
      stopCluster(cl)

    }

    var.imp <- var.imp + unlist(var.imp.list)

}

  class(var.imp) <- "varimp.CLF"

return(var.imp)
}


varimp.help.oob <- function(x, s, data, response, exposure, linear.offset, off.model){
 # browser()

  ## get strata variable for this particular tree
  sample.strata <- x$sample.strata

  ## if model contains linear offset, create linear offset for particular (possibly permuted) data
  if(linear.offset){
    ## define variable names for matrix Z
    names.Z    <- names(data)[!(names(data)%in%c(response, exposure,s))]

    ## define formula, depending on whether exposure is specified or not
    if(!is.null(exposure)){
      form.linear <- as.formula(paste0(response," ~ ", exposure, "+", paste0(names.Z,collapse="+")))
    }else{
      form.linear <- as.formula(paste0(response," ~ ", paste0(names.Z,collapse="+")))
    }

    ## define matrix Z
    Z <- model.matrix(form.linear, data = data)[,-1,drop = FALSE]

    ## catch possible errors
    if(any(!names(off.model$coefficients)== colnames(Z))){
      stop("Coefficient names of linear offset model and names of internal
           design matrix for linear offset do not coincide! Make sure, the contrast
           settings are equal to the contrast settings you used when fitting the random forest!")
    }

    ## compute offset, based on current (possibly permuted) data
    if(is.null(exposure)){
      linear.off <- Z%*%off.model$penalized
    }else{
      linear.off <- Z[,-1,drop = FALSE]%*%off.model$penalized
    }
  }else{
    linear.off <- rep(0, nrow(data))
  }

  ## extract strata numbers which are oob for this particular tree
  oob.strata <- unique(data[, names(data) == s])[which(!(unique(data[, names(data) == s]) %in% sample.strata))]

  #inbag.index <- which(data[, names(data) == s] %in% sample.strata)
  ## index of observations which are oob for this tree
  oob.index <- which(!(data[, names(data) == s] %in% sample.strata))
# browser()



  ## predict linear predictor, only for oob observations
  preds.oob <- predict(x$tree.b, newdata = data[oob.index,], response = response,
                       exposure = exposure, s = s, offset = linear.off[oob.index])


  ## split oob predictions according to strata
  preds.oob <- split(preds.oob, data[oob.index, names(data) == s])

  ## split responses of oob observations according to strata
  y.oob <- data[oob.index, names(data) == response]
  y.oob <- split(y.oob, data[oob.index, names(data) == s])


  ## compute vector of conditional likelihoods of oob strata
  oob.ll <- c()
  for(zz in seq_along(oob.strata)){
    oob.ll <- c(oob.ll, prod(exp(preds.oob[[zz]]*y.oob[[zz]]))/sum(exp(preds.oob[[zz]])))
  }


  ## return average oob conditional likelihood for this tree
  return(mean(oob.ll))

}



# varimp.help.ib <- function(x, s, data, response, exposure){
#
#   preds.oob <- predict(x$tree.b, newdata = data, response = response,
#                        exposure = exposure, s = s)
#
#
#   preds.oob <- split(preds.oob, data[oob.index, names(data) == s])
#   y.oob <- data[oob.index, names(data) == response]
#   y.oob <- split(y.oob, data[oob.index, names(data) == s])
#
#
#   oob.ll <- c()
#   for(zz in seq_along(oob.strata)){
#     oob.ll <- c(oob.ll, prod(exp(preds.oob[[zz]]*y.oob[[zz]]))/sum(exp(preds.oob[[zz]])))
#   }
#
#   return(mean(oob.ll))
#
# }


singlevar.oob <- function(var.num, model, var.names, single.trees,  n.random,
                          linear.offset, off.model){
  # browser()
  ## setze Ausgangswert für diese Variable auf 0
  var.imp <- 0

  ## wiederhole n.random mal:
  for(ooo in 1: n.random){

  ## randomisiere diese Variable
  data.rnd <- model$data
  data.rnd[, names(data.rnd) == var.names[var.num]] <- sample(data.rnd[, names(data.rnd) == var.names[var.num]])

  ## berechne Model fit pro Baum mit permutierter Variable
  var.num.list <- lapply(model$tree_list, varimp.help.oob, s = model$s, data = data.rnd,
                      response = model$response, exposure = model$exposure,
                      linear.offset = linear.offset, off.model = off.model)

  ## berechne differenz zwischen unpermutiertem fit und permutierten fit pro baum
  diff.ooo <- unlist(single.trees)-unlist(var.num.list)
  ## rechne Mittelwert dieser Differenz zu variable importance dieser variable hinzu
  var.imp <- var.imp + mean(diff.ooo)
  }

  ## teile durch Anzahl an Wiederholungen
  var.imp <- var.imp/n.random

  return(var.imp)
}


singlevar.ib <- function(var.num, model, var.names, preds.orig, strata, y, n.random){

  var.imp <- 0
  for(ooo in 1:n.random){
    data.rnd <- model$data
    data.rnd[, names(data.rnd) == var.names[var.num]] <- sample(data.rnd[, names(data.rnd) == var.names[var.num]])
    preds.ib <- predict.CLogitForest(model, newdata = data.rnd, response = model$response,
                                     exposure = model$exposure, s = model$s)

    preds.ib <- split(preds.ib, strata)

    ll.ib <- ll.orig <- c()

    for(tt in seq_along(unique(strata))){
      ll.ib <- c(ll.ib, prod(exp(preds.ib[[tt]]*y[[tt]]))/sum(exp(preds.ib[[tt]])))
      ll.orig <- c(ll.orig, prod(exp(preds.orig[[tt]]*y[[tt]]))/sum(exp(preds.orig[[tt]])))
    }

    var.imp <- var.imp + mean(ll.orig-ll.ib)/n.random
  }
    return(var.imp)

}
