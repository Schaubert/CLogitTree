

penalized.clr <- function(response,
                          stratum,
                          penalized,
                          unpenalized = NULL,
                          lambda = NULL,
                          alpha = 1,
                          p = NULL,
                          standardize = FALSE,
                          event) {

  if (missing(event) && is.factor(response)) event <- levels(response)[1]

  if (is.factor(response)) response <- (response == event) * 1

  if (!is.null(unpenalized) && !is.numeric(dim(unpenalized))) {
    unpenalized <- as.matrix(unpenalized, nrow = nrow(penalized))
  }
  if (length(p) > 0 && sum(p) != ncol(penalized)) stop("elements of p must sum to the number of penalized covariates.")



  if(is.null(lambda)){
    lambda <- find.default.lambda(response,
                                  stratum,
                                  penalized,
                                  unpenalized,
                                  alpha,
                                  p,
                                  standardize,
                                  event)
    if (is.numeric(lambda)) {lambda <- lambda[1]} else{
      lambda <- unlist(lapply(lambda, function(x) x[1]))
    }
    }else{
    if (length(lambda) > 1 && (missing(p) | is.null(p))) stop("multiple penalties are supplied in lambda, but p is missing.")
    if (length(p) != length(lambda) && !missing(p) && !is.null(p)) stop("lambda and p are not of the same length.")
    if (sum(lambda < 0) > 0) stop("lambda must be non-negative.")
  }


  if (missing(p) | is.null(p)) {
    lambda1 <- lambda
  } else {
    lambda1 <- rep(lambda, p)
  }

  if (alpha <= 0 | alpha > 1) stop("alpha must be between zero and one.")
  lambda2 <- lambda1 * (1 - alpha) / (2 * alpha)


  Y <- survival::Surv(rep(1, length(response)),
    event = (response == 1)
  )



  if (is.null(unpenalized)) {
    fit <- penalized::penalized(Y ~ strata(stratum),
      penalized = penalized,
      lambda1 = lambda1,
      lambda2 = lambda2,
      standardize = standardize
    )
  }
  else {
    fit <- penalized::penalized(Y ~ strata(stratum) + unpenalized,
      penalized = penalized,
      lambda1 = lambda1,
      lambda2 = lambda2,
      standardize = standardize
    )
  }
  return(list(
    penalized = fit@penalized,
    unpenalized = fit@unpenalized,
    converged = fit@converged,
    lambda = unique(lambda1),
    alpha = alpha,
     Y=Y))
}



subsample.clr <- function(response,
                          stratum,
                          penalized,
                          unpenalized = NULL,
                          lambda,
                          alpha = 1,
                          B = 100,
                          matB = NULL,
                          return.matB = FALSE,
                          parallel = TRUE,
                          standardize = FALSE) {
  ind.pair <- unique(stratum)
  b <- length(ind.pair)
  subsample.size <- ceiling(b / 2)

  if (is.null(matB)) {
    matB <- matrix(0, nrow = 2 * B, ncol = subsample.size)
    for (i in 1:B) {
      matB[i, ] <- sample(ind.pair,
                          size = subsample.size,
                          replace = FALSE
      )
      matB[B + i, 1:(b - subsample.size)] <- setdiff(ind.pair, matB[i, ])
    }
  }

  if (is.null(unpenalized)) {
    selB <- matrix(0, ncol = ncol(penalized), nrow = 2 * B)
  } else {
    if (!is.numeric(dim(unpenalized))) {
      unpenalized <- as.matrix(unpenalized, nrow = nrow(penalized))
    }
    selB <- matrix(0, ncol = ncol(penalized), nrow = 2 * B)
  }

  if (parallel) {
    cl <- parallel::makeCluster(getOption("cl.cores", 2), setup_timeout = 0.5)
    parallel::clusterExport(cl, varlist = c("penalized.clr"))
    selB <- t(parallel::parSapply(cl,
                                  X = 1:(2 * B),
                                  FUN = function(x,
                                                 response,
                                                 stratum,
                                                 penalized,
                                                 unpenalized,
                                                 lambda, alpha, matB,
                                                 standardize) {
                                    #require(penalized)
                                    ind <- stratum %in% matB[x, ]
                                    (penalized.clr(response[ind],
                                                   stratum[ind],
                                                   penalized = penalized[ind, ],
                                                   unpenalized = unpenalized[ind, ],
                                                   lambda = lambda,
                                                   alpha = alpha,
                                                   standardize = standardize
                                    )$penalized != 0) * 1
                                  },
                                  response, stratum, penalized,
                                  unpenalized, lambda, alpha, matB,
                                  standardize
    ))
    parallel::stopCluster(cl)
  } else {
    for (i in 1:(2 * B)) {
      ind <- stratum %in% matB[i, ]
      fit <- penalized.clr(response[ind],
                           stratum[ind],
                           penalized = penalized[ind, ],
                           unpenalized = unpenalized[ind, ],
                           lambda = lambda,
                           alpha = alpha,
                           standardize = standardize
      )
      selB[i, ] <- (fit$penalized != 0) * 1
      colnames(selB) <- names(fit$penalized)
    }
  }
  Pilambda <- colMeans(selB)
  if (return.matB) {
    return(list("Pilambda" = Pilambda, "matB" = matB))
  } else {
    return("Pilambda" = Pilambda)
  }
}
subsample.clr.v <- Vectorize(subsample.clr, vectorize.args = "lambda")






stable.clr <- function(response,
                       stratum,
                       penalized,
                       unpenalized = NULL,
                       lambda.seq = NULL,
                       alpha = 1,
                       B = 100,
                       parallel = TRUE,
                       standardize = FALSE,
                       event) {
  if (missing(event) && is.factor(response)) event <- levels(response)[1]

  if (is.factor(response)) response <- (response == event) * 1

  if (!is.null(unpenalized) && !is.numeric(dim(unpenalized))) {
    unpenalized <- as.matrix(unpenalized, nrow = nrow(penalized))
  }

  if (is.null(lambda.seq)){
    lambda.seq <- find.default.lambda(response,
                                      stratum,
                                      penalized,
                                      unpenalized,
                                      alpha,
                                      p = NULL,
                                      standardize,
                                      event)
  }

  fit <- subsample.clr(
    response = response,
    stratum = stratum,
    penalized = penalized,
    unpenalized = unpenalized,
    lambda = lambda.seq[1],
    alpha = alpha,
    B = B,
    matB = NULL,
    return.matB = TRUE,
    parallel = parallel,
    standardize = standardize
  )
  matB <- fit$matB
  if (parallel) {
    cl <- parallel::makeCluster(getOption("cl.cores", 2), setup_timeout = 0.5)
    parallel::clusterExport(cl, varlist = c("penalized.clr"))

    P <- expand.grid("B" = 1:nrow(matB), "lambda" = lambda.seq[-1])

    res <- t(parallel::parApply(cl,
                                P,
                                1,
                                FUN = function(x,
                                               response,
                                               stratum,
                                               penalized,
                                               unpenalized,
                                               matB,
                                               alpha,
                                               standardize) {
                                  #require("penalized")
                                  ind <- stratum %in% matB[x[1], ]

                                  (penalized.clr(
                                    response = response[ind],
                                    stratum = stratum[ind],
                                    penalized = penalized[ind, ],
                                    unpenalized = unpenalized[ind, ],
                                    lambda = x[2],
                                    alpha = alpha,
                                    standardize = standardize
                                  )$penalized != 0) * 1
                                },
                                response,
                                stratum,
                                penalized,
                                unpenalized,
                                matB,
                                alpha,
                                standardize
    ))
    res1 <- as.data.frame(cbind(P, res))
    res2 <- stats::aggregate(
      res,
      list(lambda = res1$lambda),
      mean
    )
    res <- t(rbind(fit$Pilambda, res2[, -c(1)]))
    parallel::stopCluster(cl)
  } else {
    res <- subsample.clr.v(
      response = response,
      stratum = stratum,
      penalized = penalized,
      unpenalized = unpenalized,
      lambda = lambda.seq[-1],
      alpha = alpha,
      B = B,
      matB = fit$matB,
      parallel = FALSE,
      standardize = standardize
    )
    res <- cbind(fit$Pilambda, res)
  }

  Pistab <- apply(res, 1, max)
  names(Pistab) <- names(fit$Pilambda)
  return(list(Pistab = Pistab, lambda.seq = lambda.seq))
}







stable.clr.g <- function(response,
                         stratum,
                         penalized,
                         unpenalized = NULL,
                         p = NULL,
                         lambda.list = NULL,
                         alpha = 1,
                         B = 100,
                         parallel = TRUE,
                         standardize = FALSE,
                         event) {
  if (missing(event) && is.factor(response)) event <- levels(response)[1]
  if (is.factor(response)) response <- (response == event) * 1
  if (!is.null(unpenalized) && !is.numeric(dim(unpenalized))) {
    unpenalized <- as.matrix(unpenalized, nrow = nrow(penalized))
  }

  if (missing(p)| is.null(p) | length(p) == 1){
    warning("valid p is not provided:
            all covariates are penalized equally.")
    if (is.null(lambda.list)) {

      temp <- stable.clr(response, stratum, penalized,
                         unpenalized, lambda.seq = NULL,
                         alpha, B, parallel, standardize, event)}else{
                           temp <- stable.clr(response, stratum, penalized,
                                              unpenalized, lambda.seq = unlist(lambda.list),
                                              alpha, B, parallel, standardize, event)
                         }
    Pistab <- temp$Pistab
    lambda.list <- temp$lambda.seq
  }else{

    g <- length(p) # the number of groups

    ind.pair <- unique(stratum)
    b <- length(ind.pair)
    subsample.size <- ceiling(b / 2)

    matB <- matrix(0, nrow = 2 * B, ncol = subsample.size)
    for (i in 1:B) {
      matB[i, ] <- sample(ind.pair,
                          size = subsample.size,
                          replace = FALSE
      )
      matB[B + i, 1:(b - subsample.size)] <- setdiff(ind.pair, matB[i, ])
    }

    if (is.null(lambda.list) | missing(lambda.list)) {
      lambda.list <- find.default.lambda(response,
                                         stratum,
                                         penalized,
                                         unpenalized,
                                         alpha,
                                         p,
                                         standardize)
    }
    if (parallel) {
      cl <- parallel::makeCluster(getOption("cl.cores", 2), setup_timeout = 0.5)
      parallel::clusterExport(cl, varlist = c("penalized.clr"))

      subslist <- c(list(B = 1:(2 * B)), lambda.list)
      P <- expand.grid(subslist)

      res <- t(parallel::parApply(cl,
                                  P,
                                  1,
                                  FUN = function(x,
                                                 response, stratum, penalized, unpenalized,
                                                 p, matB, alpha, standardize) {
                                    ind <- stratum %in% matB[x[1], ]
                                    (penalized.clr(
                                      response = response[ind],
                                      stratum = stratum[ind],
                                      penalized = penalized[ind, ],
                                      unpenalized = unpenalized[ind, ],
                                      p = p,
                                      lambda = as.numeric(x[-1]),
                                      alpha = alpha,
                                      standardize = standardize)$penalized != 0) * 1
                                  },
                                  response,
                                  stratum,
                                  penalized,
                                  unpenalized,
                                  p,
                                  matB,
                                  alpha,
                                  standardize))

      res1 <- as.data.frame(cbind(P, res))
      res2 <- stats::aggregate(res, by = as.list(res1[, 2:(g + 1)]), mean)
      res <- t(res2[, -c(1:g)])
      parallel::stopCluster(cl)
    } else {
      for (i in 1:nrow(P)) {
        ind <- stratum %in% matB[P[i, 1], ]
        fit <- penalized.clr(response[ind],
                             stratum[ind],
                             penalized = penalized[ind, ],
                             unpenalized = unpenalized[ind, ],
                             lambda = as.numeric(P[i, -1]),
                             alpha = alpha,
                             standardize = standardize )
        selB[i, ] <- (fit$penalized != 0) * 1
        colnames(selB) <- names(fit$penalized)
      }
    }

    Pistab <- apply(res, 1, max)
    names(Pistab) <- colnames(penalized)}

  return(list(Pistab = Pistab, lambda.list = lambda.list))
}









find.default.lambda <- function(response, stratum, penalized,
                                unpenalized = NULL,
                                alpha = 1,
                                p = NULL,
                                standardize = FALSE,
                                event,
                                nfolds = 10){

  if (missing(event) && is.factor(response)) event <- levels(response)[1]

  if (is.factor(response)) response <- (response == event) * 1

  if (!is.null(unpenalized)) {X <- cbind(unpenalized, penalized)
  nc <- ncol(unpenalized)} else{
    X <- penalized
    nc <- 0
  }

  if(standardize == T) X <- apply(X, 2, function(x)
    x/sqrt((length(x)-1)/length(x)*var(x)))

  if(missing(p) | is.null(p)){
    lambda.seq <- default.lambda(X, response, stratum, alpha, nfolds = nfolds)
  } else{
    g <- length(p)
    le <- c(1, p[-g]+1) + nc
    ue <- cumsum(p) + nc
    lambda.seq <- list()
    for(i in 1:g){
      lambda.seq[[i]] <- default.lambda(X = X[, c((nc != 0)*(1:nc),le[i]:ue[i])],
                                        response,
                                        stratum,
                                        alpha,
                                        nfolds = nfolds)
    }
  }

  return("lambda.seq" = lambda.seq)
}






default.lambda <- function(X, Y, stratum, nfolds = 10, alpha = 1){
  fit <- clogitL1::clogitL1(X, Y, stratum, alpha = alpha)
  cvfit <- clogitL1::cv.clogitL1(fit, numFolds = nfolds)
  d_lambda <- unique(c(
    exp((3*cvfit$minCV_lambda - cvfit$minCV1se_lambda)/2),
    exp(cvfit$minCV_lambda),
    exp((cvfit$minCV_lambda + cvfit$minCV1se_lambda)/2),
    exp(cvfit$minCV1se_lambda)))
  return(d_lambda)
}


