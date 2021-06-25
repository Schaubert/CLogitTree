prune.CLogitTree <- function(x,
                             alpha){


  alpha_adj <- alpha/ncol(x$Z)


  if(!any(x$pvalue>alpha_adj)){
    iter <- length(x$pvalue)
  }else{
    iter <- min(which(x$pvalue>alpha_adj))
  }


  model      <- x$model[[iter]]
  params     <- x$param[[iter]]
  params_fit <- x$param_fit[[iter]]
  design     <- x$design[[iter]]
  pvalues    <- x$pvalue[1:iter]
  devs       <- x$dev[1:iter]
  crits      <- x$crit[1:iter]

  if(iter==1){
    if(length(x$exposure)>0){
      beta_hat <- unlist(model$penalized)
    } else{
      beta_hat <- NULL
    }
    gamma_hat <- NULL
    splits   <- NULL
  }

  if(iter>1){
    if(length(x$exposure)>0){
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
    splits <- x$splits[1:(iter-1),]
  }

  to_return <-  list("beta_hat"=beta_hat,
                     "gamma_hat"=gamma_hat,
                     "splits"=splits,
                     "Z"=x$Z,
                     "y"=x$y,
                     "model"=model,
                     "design"=design,
                     "param"=params,
                     "pvalue"=pvalues,
                     "dev"=devs,
                     "crit"=crits,
                     "call"=x$call)

  class(to_return) <- "CLogitTree"
  return(to_return)
}
