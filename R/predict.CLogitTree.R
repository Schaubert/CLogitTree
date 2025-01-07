#' Predict function for CLogitTree
#'
#' Predicts the linear predictor of estimated CLogitTree models
#'
#' @param object CLogitTree object
#' @param response Name of response variable (character)
#' @param exposure Name of (optional) separate exposure variable (character)
#' @param s Name of strata variable (character)
#' @param newdata (Optional) new data frame. If \code{newdata = NULL}, the original data set from \code{object} is used.
#' @param type Type of prediction that is required. For option \code{"link"} (default), the prediction of the
#' linear predictor per observation is returned. For option \code{"response"} (default), the prediction of the
#' probability per observation is returned. This probability is conditional on the stratum,
#' i.e. the probabilities per stratum add up to 1. For option \code{"loglik"}, the prediction of the conditional log-likelihood
#' per stratum is returned. For option \code{"tree"}, the prediction of the tree (without exposure effect)
#' per observation is returned.
#' @param offset Optional offset argument, mainly for internal use
#' @param ... Further predict arguments
#' @return Vector containing predicted linear predictors
#' @author Gunther Schauberger: \email{gunther.schauberger@@tum.de} \cr
#' Moritz Berger: \email{moritz.berger@@imbie.uni-bonn.de}
#' @seealso \code{\link{CLogitTree}}
#' @export
predict.CLogitTree <- function(object,
                               response,
                               exposure=NULL,
                               s,
                               newdata=NULL,
                               type = c("link", "response", "loglik", "tree"),
                               offset = NULL,
                               ...){
  type <- match.arg(type)


 if(is.null(exposure)){
   if(!is.character(object$exposure.name)){
      exposure <- NULL
   }else{
     exposure <- object$exposure.name
   }
 }

  if(is.null(newdata)){

    Z <- object$Z
    nvar <- ncol(Z)

      newdata <- object$data
      strata <- object$strata
      y <- object$y


  } else{
    strata <- newdata[, names(newdata) == s]


    y <- newdata[, names(newdata) == response]
    Z <- newdata[,!names(newdata)%in%c(response,exposure,s)]
    nvar <- ncol(Z)

    if(!is.null(names(Z))){
      var_names <- names(Z)
    } else{
      var_names <- paste0("x",1:nvar)
    }

    for(j in 1:nvar){
      if(is.factor(Z[,j])){
        if(!is.ordered(Z[,j])){
          Z[,j] <- mod_factors(y,Z[,j])
        } else{
          Z[,j] <- as.integer(Z[,j])
        }
      }
    }

  }

  gamma_hat <- object$gamma_hat_sym

  if(!is.null(gamma_hat)){
    ordered_values <- lapply(1:nvar, function(j) ord_values(object$Z[,j]))
    n_levels       <- sapply(ordered_values,length)
    thresholds     <- lapply(ordered_values,thresh)
    n_s            <- n_levels-1
    design_lower <- designlists(Z,nvar,n_s,n_levels,ordered_values)[[1]]
    design_lower <- do.call(cbind, design_lower)
    design_upper <- designlists(Z,nvar,n_s,n_levels,ordered_values)[[2]]
    design_upper <- do.call(cbind, design_upper)

    DM <- as.data.frame(cbind(design_lower, design_upper, "int"=1))
    form <- formula(paste("~",paste(names(gamma_hat), collapse="+"),"-1"))
    mm <- model.matrix(form, data=DM)
    colnames(mm) <- check_names(colnames(mm), names(gamma_hat))

    eta_gamma <- mm%*%gamma_hat[colnames(mm)]
  } else{
    eta_gamma <- rep(0,nrow(Z))
  }


  if((!is.null(exposure)) & (type != "tree")){
    beta_hat  <- object$beta_hat
    if(is.null(newdata)){
      expo <- object$exposure
    } else{
      expo <- newdata[, names(newdata) == exposure]
    }
    expo <- model.matrix(~expo)[,-1,drop=FALSE]
    eta <- eta_gamma + expo%*%beta_hat
  } else{
    eta <- eta_gamma
  }


  if(!is.null(offset)){
    if(length(offset)==1){
      rep(offset, length(eta))
      warning("Offset value is scalar, but should have the length of data or newdata.")
    }

    if(length(offset)==length(eta)){
      eta <- eta + offset
    }else{
      stop("Length of offset does not coincide with length of eta!")
    }

  }


  if(type == "loglik"){


    pred.list <- split(eta, strata)
    y.list <- split(y, strata)


    ll <- c()
    for(zz in 1:length(y.list)){
      ll <- c(ll, prod(exp(pred.list[[zz]]*y.list[[zz]]))/sum(exp(pred.list[[zz]])))
    }
    return(ll)
  }else{
    if(type == "response"){
      pred.list <- split(eta, strata)
      y.list <- split(y, strata)


      mu <- c()
      for(zz in 1:length(y.list)){
        mu <- c(mu, exp(pred.list[[zz]])/sum(exp(pred.list[[zz]])))
      }
      return(mu)
    }else{
      return(c(eta))
    }
  }



}
