#' Predict function for CLogitForest
#'
#' Predicts the linear predictor of estimated CLogitForest models
#'
#' @param object CLogitForest object
#' @param response Name of response variable (character)
#' @param exposure Name of (optional) separate exposure variable (character)
#' @param s Name of strata variable (character)
#' @param newdata (Optional) new data frame. If \code{newdata = NULL}, the original data set from \code{object} is used.
#' @param type Type of prediction that is required. For option \code{"link"} (default), the prediction of the
#' linear predictor per observation is returned. For option \code{"response"} (default), the prediction of the
#' probability per observation is returned. This probability is conditional on the stratum,
#' i.e. the probabilities per stratum add up to 1. For option \code{"loglik"}, the prediction of the conditional log-likelihood
#' per stratum is returned. For option \code{"forest"}, the prediction of the random forest (without exposure effect)
#' per observation is returned.
#' @param ncores Number of parallel cores to use
#' @param median.exposure Shall the median of the exposure effects from single trees be used as
#' estimate (instead of arithmetic mean)
#' @param ... Further predict arguments
#' @return Vector containing predicted linear predictors
#' @author Gunther Schauberger: \email{gunther.schauberger@@tum.de} \cr
#' Moritz Berger: \email{moritz.berger@@imbie.uni-bonn.de}
#' @seealso \code{\link{CLogitForest}}
#' @export
predict.CLogitForest <- function(object,
                               response = NULL,
                               exposure = NULL,
                               s = NULL,
                               newdata = NULL,
                               ncores = 1,
                               type = c("link", "response", "loglik", "forest"),
                               median.exposure = FALSE,
                               ...){

  type <- match.arg(type)

  linear.offset <- object$linear.offset


  if(is.null(newdata)){
    newdata <- object$data
    response <- object$response
    exposure <- object$exposure
    s <- object$s
  }
 # browser()

  if(linear.offset){
    names.Z    <- names(newdata)[!(names(newdata)%in%c(response, exposure,s))]
    # names.Z <- names(object$lin.off.model$coefficients)
    # if(!is.null(exposure)){
    #   names.Z <- names.Z[-1]
    # }
    if(!is.null(exposure)){
      form.linear <- as.formula(paste0(response," ~ ", exposure, "+", paste0(names.Z,collapse="+")))
    }else{
      form.linear <- as.formula(paste0(response," ~ ", paste0(names.Z,collapse="+")))
    }

    Z <- model.matrix(form.linear, data = newdata)[,-1,drop = FALSE]

    # Z    <- newdata[,!(names(newdata)%in%c(response, exposure,s))]

    if(any(!names(object$lin.off.model$coefficients)== colnames(Z))){
      stop("Coefficient names of linear offset model and names of internal
           design matrix for linear offset do not coincide! Make sure, the contrast
           settings are equal to the contrast settings you used when fitting the random forest!")
    }

    if(is.null(exposure)){
      linear.off <- Z%*%object$lin.off.model$penalized
    }else{
      linear.off <- Z[,-1,drop = FALSE]%*%object$lin.off.model$penalized
    }
  }else{
    linear.off <- rep(0, nrow(newdata))
  }


  if(ncores == 1){
    pred <- lapply(object$tree_list, help.pred,
                        newdata = newdata,
                        s = s,
                        exposure = exposure,
                        response = response,
                   #offset = offset,
                   linear.off = linear.off)
  }else{
    cl <- makeCluster(ncores, outfile = "")

    clusterExport(cl, varlist = c("newdata", "response", "exposure", "s", "linear.off"),
                  envir = sys.frame(sys.nframe()))

    pred <- parLapply(cl, object$tree_list, help.pred,
                                 newdata = newdata,
                                 s = s,
                                 exposure = exposure,
                                 response = response,
                      #offset = offset,
                      linear.off = linear.off)
    stopCluster(cl)
  }


  pred.mat <- matrix(unlist(pred),ncol = object$ntree)

# browser()

    pred <- rowMeans(pred.mat)

    if(type != "forest"){
      if(!is.null(exposure)){
      if(median.exposure){
        pred <- pred + object$beta_median_hat * newdata[,exposure]
      }else{
        pred <- pred + object$beta_hat * newdata[,exposure]
      }
      }
    }


  # if(offset){
  #   if(type!="forest"){
  #
  #     newdata$offset.vec <- pred
  #     offset.model <- clogistic(as.formula(paste0(response," ~ ", exposure, "+ offset(offset.vec)")),
  #                            strata = newdata[,s], data = newdata)
  #
  #   pred <- newdata[,exposure]*coef(offset.model)[1] + newdata$offset.vec
  #   }
  # }

    # browser()

  if(type == "loglik"){
    pred.list <- split(pred, newdata[, names(newdata) == s])
    y <- newdata[, names(newdata) == response]
    y.list <- split(y, newdata[, names(newdata) == s])


    ll <- c()
    for(zz in 1:length(y.list)){
      ll <- c(ll, prod(exp(pred.list[[zz]]*y.list[[zz]]))/sum(exp(pred.list[[zz]])))
    }
    return(ll)
  }else{
    if(type == "response"){
      pred.list <- split(pred, newdata[, names(newdata) == s])
      y <- newdata[, names(newdata) == response]
      y.list <- split(y, newdata[, names(newdata) == s])


      mu <- c()
      for(zz in 1:length(y.list)){
        mu <- c(mu, exp(pred.list[[zz]])/sum(exp(pred.list[[zz]])))
      }
      return(mu)
    }else{
      return(pred)
      }
  }
}

### old version of help.pred with old version of offset
#   help.pred <- function(x, newdata, s, exposure, response, offset, type){
#     if(offset){
#       return <- predict(x$tree.b, newdata = newdata,
#               s = s, response = response, exposure = exposure, type = "tree")
#       }else{
#         type.x <- "link"
#         if(type == "forest"){
#           type.x <- "tree"
#         }
#         return <- predict(x$tree.b, newdata = newdata,
#             s = s, response = response, exposure = exposure, type = type.x)
#       }
#
#     if(length(return)==1){return <- rep(0,nrow(newdata))}
#
#     return
# }

### new help.pred version
help.pred <- function(x, newdata, s, exposure, response, linear.off){


    return <- predict(x$tree.b, newdata = newdata,
                      s = s, response = response, exposure = exposure, type = "tree", offset = linear.off)


  if(length(return)==1){return <- rep(0,nrow(newdata))}

  return
}
