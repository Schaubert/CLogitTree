prune <- function (x, ...) {
  UseMethod("prune", x)
}

prune.CLogitTree <- function(x,
                             alpha){

  alpha_adj <- alpha/ncol(x$Z)
  iter      <- which.min(x$pvalue>alpha_adj)

  model  <- x$model[[iter]]
  params <- x$param[[iter]]
  design <- x$design[[iter]]
  pvalue <- x$pvalue[[iter]]
  dev    <- x$devs[[iter]]
  crit   <- x$crits[[iter]]

  beta_hat <- coefficients(model)[1]

  if(iter>1){
    gamma_hat <- coefficients(model)[-1]
    gamma_hat[is.na(gamma_hat)] <- 0
    gamma_hat <- gamma_hat[params]
    splits <- x$splits[1:(iter-1),]
  } else{
    gamma_hat <- NULL
    splits   <- NULL
  }

  to_return <-  list("beta_hat"=beta_hat,
                     "gamma_hat"=gamma_hat,
                     "splits"=splits,
                     "Z"=x$Z,
                     "y"=x$y,
                     "model"=model,
                     "design"=design,
                     "param"=params,
                     "pvalue"=pvalue,
                     "dev"=dev,
                     "crit"=crit,
                     "call"=x$call)

  class(to_return) <- "CLogitTree"
  return(to_return)
}
