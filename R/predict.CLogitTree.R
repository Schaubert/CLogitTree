#' Predict for CLogitTree
#' @export
predict.CLogitTree <- function(object,
                               response,
                               exposure=NULL,
                               s,
                               newdata=NULL,
                               ...){

  if(is.null(newdata)){

    Z <- object$Z
    nvar <- ncol(Z)

  } else{

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

  gamma_hat <- object$gamma_hat

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
    eta_gamma <- 0
  }

  if(!is.null(exposure)){
    beta_hat  <- object$beta_hat
    if(is.null(newdata)){
      expo <- object$exposure
    } else{
      expo <- newdata[, names(newdata) == exposure]
    }
    if(is.factor(expo)){
      expo <- as.numeric(expo)-1
    }
    eta <- eta_gamma + expo*beta_hat
  } else{
    eta <- eta_gamma
  }

 return(eta)
}
