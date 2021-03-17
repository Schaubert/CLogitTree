cv.CLogitTree <- function(x, K= 5, p.range = seq(0.2, 0.01, by = -0.01)){
browser()
  n.p <- length(p.range)

  fold.size <- floor(length(unique(x$strata))/K)
  size.vec <- rep(fold.size, K) +  c(rep(1, length(unique(x$strata))%%K), rep(0, K-length(unique(x$strata))%%K))
  strata.list <- split(sample(unique(x$strata)), as.factor(unlist(strsplit(paste0(strrep(1:K, size.vec), collapse = ""),""))))

  X <- cbind(x$Z, x$strata, x$exposure)
  names(X)[(ncol(X)-1):ncol(X)] <- c("sk", "expk")

  lp.vec <- rep(0, length(p.range))

  for(k in 1:K){
    index.k <- which(X$sk %in% strata.list[[k]])
    index.k.not <- which(!(X$sk %in% strata.list[[k]]))

    X.k <- X[index.k.not,]
    y.k <- y[index.k.not]

    ## falsche design matrix, muss man selber bauen
    design.k <- x$Z[index.k,]
    exposure.k <- x$exposure[index.k]
    strata.k <- x$strata[index.k]

    tree.k <- CLogitTree(y.k, X.k, "expk", "sk", alpha = max(p.range), nperm = x$nperm, minnodesize=x$minnodesize, perm_test=TRUE, trace=x$trace, fit=TRUE)

    for(pp in seq_along(p.range)){
      tree.k.pp <- prune(tree.k, p.range[pp])

      eta.k.pp <- design.k%*%tree.k.pp$gamma_hat + exposure.k*tree.k.pp$beta.hat
      eta.list <- split(eta.k.pp, strata.k)
      lp.vec[pp] <- lp.vec[pp] + prod(lapply(eta.list, function(x){exp(x[1])/sum(exp(x))}))
    }

    p.final <- p.range[which.min(lp.vec)]

    to_return <- prune(x, p.final)

  }

  class(to_return) <- "CLogitTree"
  return(to_return)
}
