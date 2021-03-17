

CLogitForest <- function(y,
                       X,
                       exposure,
                       s,
                       alpha = 0.25,
                       nperm = 100,
                       B = 20,
                       mtry = 3,
                       subsample = 0.632,
                       trace=TRUE){



  Z    <- X[,!names(X)%in%c(exposure,s)]
  n    <- length(y)
  nvar <- ncol(Z)

  strata <- unique(X[,names(X)%in%c(s)])
  n.strata <- length(strata)

exp.effect <- c()

  for(b in 1:B){

    sample.strata <- sample(strata, ceiling(subsample*n.strata), replace = FALSE)

    subset.index <- which(X[,names(X)%in%c(s)] %in% sample.strata)

    y.b <- y[subset.index]

    var.forest <- sample(1:nvar, mtry)
    X.b <- cbind(X[subset.index,names(X)%in%c(exposure,s)], Z[subset.index,var.forest])

    tree.b <- CLogitTree(y = y.b,
               X = X.b,
               exposure = exposure,
               s = s,
               alpha = alpha,
               nperm = nperm,
               trace = trace)

exp.effect <- c(exp.effect, tree.b$beta_hat)
cat("\n Random Forest",b,"out of",B,"\n")
  }

  to_return <- (list("beta_hat"=mean(exp.effect),
                     "call"=match.call()))

  class(to_return) <- "CLogitForest"
  return(to_return)

}
