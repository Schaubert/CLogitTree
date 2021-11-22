#
#
# boot.CLfun <- function(data, index, model){
#
#   X = cbind(model$Z, model$exposure, model$strata)[index,]
#   names(X)[(ncol(X)-1):ncol(X)] <- c("Xexpo", "Xstrata")
#
#   y <- model$y[index]
#
#   ret <- CLogitTree(y = y,
#              X = X,
#              exposure="Xexpo",
#              s = "Xstrata",
#              alpha = model$alpha,
#              nperm = model$nperm,
#              minnodesize=model$minnodesize,
#              perm_test=model$perm_test,
#              mtry=model$mtry,
#              lambda=model$lambda,
#              trace=model$trace,
#              fit=TRUE)
#
# ret$beta_hat
#
# }



bootci.CLogitTree <- function(model, B = 100, alpha = 0.05){

  X = cbind(model$Z, model$exposure, model$strata, model$y)
  names(X)[(ncol(X)-2):ncol(X)] <- c("Xexpo", "Xstrata","Xy")

  index <- levels(droplevels(as.factor(model$strata)))

  if(is.null(model$beta_hat)){
    beta_hat <- CLogitTree(data = X,
                       response = "Xy",
                       exposure="Xexpo",
                       s = "Xstrata",
                       alpha = model$alpha,
                       nperm = model$nperm,
                       minnodesize=model$minnodesize,
                       perm_test=model$perm_test,
                       mtry=model$mtry,
                       lambda=model$lambda,
                       trace=model$trace,
                       fit=TRUE)$beta_hat
  }else{
    beta_hat <- model$beta_hat
  }



  beta.hat.vec <- rep(NA,B)

  for(j in 1:B){
    cat("Starting bootstrap iteration ",j,"out of ",B,"\n")
  index2 <- sample(index, size = length(index), replace = TRUE)



  X2 <- c()
  for(i in 1:length(index2)){
    Xnew <- X[X$Xstrata == index2[i],]
    Xnew$Xstrata <- rep(i, nrow(Xnew))
    X2 <- rbind(X2, Xnew)
  }

  ret <- CLogitTree(data = X2,
                    response = "Xy",
                    exposure="Xexpo",
                    s = "Xstrata",
                    alpha = model$alpha,
                    nperm = model$nperm,
                    minnodesize=model$minnodesize,
                    perm_test=model$perm_test,
                    mtry=model$mtry,
                    lambda=model$lambda,
                    trace=model$trace,
                    fit=TRUE)

  beta.hat.vec[j] <- ret$beta_hat
  }


  return(list(beta_hat = beta_hat, ci_beta = quantile(beta.hat.vec, c(alpha/2, 1-alpha/2))))
}






# boot.CLogitTree(data, model, R){
#
#   boot(data, boot.CLfun, R = R, model = model)
# }
#
# set.seed(1860)
# b1 <- boot.CLogitTree2(sub.mat2, tree.mod2)
#
# set.seed(1860)
# boot::boot(sub.mat2, boot.CLfun, R = 10, model = tree.mod2)
