#' Function to perform CLogitTree
#'
#' Performs CLogitTree, a tree-based method for the analysis of matched case-control studies. The estimation is embedded into the framework of conditional logistic regression.
#' Splits are incorporated into the model as separate dummy variables specifying the terminal nodes. Additional methods for pruning, plotting and bootstrap are available.
#'
#' @param data Data frame (containing all required variable)
#' @param response Name of response variable (character)
#' @param exposure Name of (optional) separate exposure variable (character)
#' @param s Name of strata variable (character)
#' @param alpha Level of significance for internal permutation tests
#' @param nperm Number of permutations in internal permutation tests
#' @param minnodesize Minimal node size in order to be eligible for further splitting
#' @param minbucket Minimum number of observations in any terminal node
#' @param depth_max Maximum depth of the tree, with the root node counted as depth 0. If \code{NULL} (default), the size of the trees is not restricted.
#' @param perm_test \code{TRUE} for regular use where  permutation tests are done.
#' @param mtry Solely for internal use
#' @param lambda Tuning parameter for optional L2 penalty
#' @param offset Optional offset used for model fitting
#' @param print.trace Shall trace of permutation tests be printed?
#' @param ncores Number of parallel nodes to use for permutation tests.
#' @param fit Shall the internally fitted models be returned (required for the use of \code{\link[CLogitTree]{prune}})?
#' @param epsilon Convergence tolerance. Iteration continues until the relative change
#' in the conditional log likelihood is less than eps. Must be positive.
#'
#' @return
#' \item{beta_hat}{Estimate for separate exposure  effect}
#' \item{gamma_hat}{Estimates for terminal nodes  (with reference node)}
#' \item{gamma_hat_sym}{Estimates for terminal nodes (with symmetric side constraint)}
#' \item{splits}{Detailed information about all executed splits during the fitting process (matrix)}
#' \item{Z}{Design matrix of all confounding variables}
#' \item{y}{Outcome variable}
#' \item{y_tab}{Table of the number of cases and the number of controls in each node}
#' \item{strata}{Strata variable}
#' \item{exposure}{Exposure variable}
#' \item{model}{Models fitted internally}
#' \item{design}{Internally built design matrix including all indicator variables}
#' \item{param}{Fitted coefficients (solely for internal use)}
#' \item{param_fit}{Fitted coefficients (solely for internal use)}
#' \item{pvalue}{P-values of the permutation tests performed during the fitting process}
#' \item{dev}{Maximal value statistics of the selected variables in each iteration during the fitting process}
#' \item{crit}{Critical values of the permutation tests during the fitting process}
#' \item{BIC}{BICs of the internally fitted models (based on the unpenalized likelihood)}
#' \item{minnodesize}{Minimal node size in order to be eligible for further splitting}
#' \item{minbucket}{Minimum number of observations in any terminal node}
#' \item{depth_max}{Maximum depth of the tree}
#' \item{lambda}{L2 tuning parameter}
#' \item{print.trace}{print.trace argument from call}
#' \item{nperm}{nperm argument from call}
#' \item{alpha}{Level of significance for permutation tests}
#' \item{perm_test}{perm_test argument from call}
#' \item{mtry}{mtry argument from call}
#' \item{ncores}{Number of parallel nodes to use for permutation tests.}
#' \item{call}{function call}
#' \item{prunedBIC}{For internal use}
#' \item{epsilon}{Convergence tolerance}
#' @author Gunther Schauberger: \email{gunther.schauberger@@tum.de} \cr
#' Moritz Berger: \email{Moritz.Berger@@imbie.uni-bonn.de}
#' @seealso \code{\link{plot.CLogitTree}}, \code{\link{bootci.CLogitTree}}, \code{\link{prune.CLogitTree}}, \code{\link{predict.CLogitTree}}
#' @examples
#' data(illu.small)
#'
#' set.seed(1860)
#' illu.tree <- CLogitTree(illu.small, response = "y", exposure = "x", s = "strata",
#'                         alpha = 0.05, nperm = 20)
#'
#' plot(illu.tree)
#' ## preferred procedure: pruning using BIC
#' set.seed(1860)
#' illu.tree <- CLogitTree(illu.small, response = "y", exposure = "x", s = "strata",
#'                         perm_test = FALSE, depth_max=4)
#' illu.tree <- pruneBIC(illu.tree)
#' @export
CLogitTree <- function(data,
                       response,
                       exposure = NULL,
                       s,
                       alpha = 0.05,
                       nperm = 500,
                       minnodesize = 10,
                       minbucket = 5,
                       depth_max = NULL,
                       perm_test = TRUE,
                       mtry = NULL,
                       lambda = 0,
                       offset = NULL,
                       print.trace = FALSE,
                       ncores = 1,
                       fit = TRUE,
                       epsilon = 1e-3){

  y <- data[, names(data) == response]
  Z    <- data[,!names(data)%in%c(response, exposure,s)]

  # exclude variables from Z with only one value
  exc <- numeric(ncol(Z))
  for(j in 1:ncol(Z)){
    if(length(unique(Z[,j]))==1){
      exc[j] <- 1
    }
  }
  Z <- Z[ ,exc==0]

  n    <- length(y)
  nvar <- ncol(Z)

    #browser()

  if(!is.null(names(Z))){
    var_names <- names(Z)
  } else{
    var_names <- paste0("x",1:nvar)
  }

  # offset
  if(is.null(offset)){
    offset_vec <- rep(0, n)
  }else{
    if(length(offset)==1){
      offset_vec <- rep(offset, n)
    } else{
      offset_vec <- offset
    }
  }

  # modify design
  for(j in 1:nvar){
    if(is.factor(Z[,j])){
      if(!is.ordered(Z[,j])){
        Z[,j] <- mod_factors(y,Z[,j])
      } else{
        Z[,j] <- as.integer(Z[,j])
      }
    }
  }

  ordered_values <- lapply(1:nvar, function(j) ord_values(Z[,j]))
  n_levels       <- sapply(ordered_values,length)
  thresholds     <- lapply(ordered_values,thresh)
  n_s            <- n_levels-1

  mod_potential <- list()
  devs          <- c()
  crits         <- c()
  pvalues       <- c()
  BICs          <- c()
  splits        <- data.frame("variable"=numeric(),
                              "split"=numeric(),
                              "level"=numeric(),
                              "node"=numeric(),
                              "threshold"=numeric(),
                              "number"=numeric(),
                              "left"=numeric(),
                              "right"=numeric(),
                              "association"=numeric(),stringsAsFactors=FALSE)

  params        <- list()
  params_fit    <- list()
  vars_evtl     <- list()
  splits_evtl   <- list()
  which_obs     <- list()
  y_tab         <- list()
  numbers       <- list()
  design        <- list()
  count         <- 1

  params[[1]]      <- "int"
  which_obs[[1]]   <- matrix(1:n,nrow=1)
  y_tab[[1]]       <- matrix(table(y), nrow=1)
  vars_evtl[[1]]   <- matrix(1:nvar,nrow=1)
  splits_evtl[[1]] <- lapply(1:nvar, function(var) matrix(1:n_s[var],nrow=1))
  numbers[[1]]     <- 1

  dat0   <- design[[1]] <- data.frame("y"=y,"int"=rep(1,n),data)
  if(!is.null(exposure)){
    form0 <- formula(paste0("y~", exposure, "+ offset(offset_vec)"))
    mod0   <- clogistic(form0, strata = dat0[,s], data = dat0, eps = epsilon, toler.chol = epsilon/10)
    BICs[1] <- (-2)*mod0$loglik[2]+log(n)*length(coef(mod0))
    mm <- model.matrix(mod0, data=dat0)
    invisible(capture.output(mod_potential[[1]] <- penalized.clr(response=dat0$y, stratum=dat0[,s], penalized=mm[,2, drop=FALSE], unpenalized=NULL,
                                                                   offset=offset_vec, lambda=0, alpha=10e-20, epsilon = epsilon)))
    params_fit[[1]] <- names(coefficients(mod0))
  } else{
    mod0   <- mod_potential[[1]] <- NULL
    params_fit[[1]] <- NULL
  }

  design_lower <- designlists(Z,nvar,n_s,n_levels,ordered_values)[[1]]
  design_upper <- designlists(Z,nvar,n_s,n_levels,ordered_values)[[2]]
  sig      <- TRUE
  anysplit <- TRUE

  # minimal node size > n?
  if(length(which(!is.na(which_obs[[1]][1,])))<minnodesize){
    for(var in 1:nvar){splits_evtl[[1]][[var]][1,] <- rep(NA, n_s[var])}
  }
  # depth_max==0?
  if(!is.null(depth_max)){
    if(depth_max==0){
      for(var in 1:nvar){splits_evtl[[1]][[var]][1,] <- rep(NA, n_s[var])}
    }
  }
  anysplit <- !all(is.na(unlist(splits_evtl[[1]])))

  while(sig & anysplit){

    # minbucket
    for(var in 1:nvar){
      n_knots <- length(params[[count]])
      for(kn in 1:n_knots){
        splits_aktuell <- splits_evtl[[count]][[var]][kn,]
        if(length(splits_aktuell)>0){
          for(j in splits_aktuell){
            thresh <- ordered_values[[var]][1:n_s[var]][j]
            obs0   <- which_obs[[count]][kn,]
            obs1   <- obs2 <- obs0
            obs1[Z[,var]>thresh] <- NA
            obs2[Z[,var]<=thresh] <- NA
            if(length(which(!is.na(obs1)))<minbucket | length(which(!is.na(obs2)))<minbucket){
              splits_evtl[[count]][[var]][kn,j] <- NA
            }
          }
        }
      }
    }
    anysplit <- !all(is.na(unlist(splits_evtl[[count]])))

    if(anysplit){
# browser()
      # draw mtry variables
      n_knots <- length(params[[count]])
      if(!is.null(mtry)){
        tryvars <- sapply(1:n_knots, function(kn) sample(1:nvar, mtry))
        tryvars <- matrix(tryvars, ncol=n_knots)
      } else{
        tryvars <- sapply(1:n_knots, function(kn) 1:nvar)
      }

      # estimate all models
      dv <- lapply(1:nvar,function(var){
                deviances <- matrix(rep(0,n_s[var]*n_knots),ncol=n_knots)
                which_knots <- which(!is.na(vars_evtl[[count]][,var]))
                for(kn in which_knots){
                  if(var %in% tryvars[,kn]){
                    deviances[,kn] <- allmodels(var,exposure,s,kn,count,design_lower,design_upper,splits_evtl,params,dat0,mod0,n_s,offset_vec,epsilon)
                  }
                }
                return(deviances)
            })

      if(max(unlist(dv))>0){

        # select optimum
        variable <- which.max(lapply(1:nvar,function(j) max(dv[[j]])))
        split    <- as.numeric(which(dv[[variable]]==max(dv[[variable]]),arr.ind=TRUE)[,1])
        knoten   <- as.numeric(which(dv[[variable]]==max(dv[[variable]]),arr.ind=TRUE)[,2])
        if(length(split)>1){
          split  <- split[1]
          knoten <- knoten[1]
        }
        param_old   <- params[[count]][knoten]
        level       <- length(strsplit(param_old,":")[[1]])
        number      <- numbers[[count]][knoten]
        left        <- max(numbers[[count]])+1
        right       <- max(numbers[[count]])+2

        param_new   <- paste(param_old,c(colnames(design_lower[[variable]])[split],colnames(design_upper[[variable]])[split]),sep=":")


        # compute permutation test
        if(perm_test){


          ## auskommentiert am 22.03.22 von GS, um parallele Implementierung zu versuchen
          # dev <- rep(NA,nperm)
          # for(perm in 1:nperm){
          #   dev[perm] <- one_permutation(variable,exposure,s,knoten,count,nvar,n_levels,ordered_values,
          #                                Z,which_obs,splits_evtl,params,dat0,mod0,n_s)
          #   if(print.trace){
          #     cat(".")
          #   }
          # }

          seeds <- abs(round(rnorm(nperm) * 1e8))

          if(ncores ==1){
            dev <- sapply(seeds, one_permutation2,
                          var = variable,
                          exposure = exposure,
                          s = s,
                          kn = knoten,
                          count = count,
                          nvar = nvar,
                          n_levels = n_levels,
                          ordered_values = ordered_values,
                          DM_kov = Z,
                          which_obs = which_obs,
                          splits_evtl = splits_evtl,
                          params = params,
                          dat0 = dat0,
                          mod0 = mod0,
                          n_s = n_s,
                          offset = offset_vec,
                          epsilon = epsilon)
          }else{
            # libspath <- .libPaths()
            # cl <- makeCluster(ncores, outfile = "CLogitTree_log.txt")
            cl <- makeCluster(ncores, outfile = "")

            clusterExport(cl, varlist = c("variable","exposure", "s", "knoten", "count",
                                          "nvar", "n_levels", "ordered_values", "Z",
                                          "which_obs", "splits_evtl", "params", "dat0",
                                          "mod0", "n_s","designlists","allmodels",
                                          "one_model","epsilon","offset_vec"),
                          envir = sys.frame(sys.nframe()))
            # clusterEvalQ(cl, .libPaths(libspath))
            clusterEvalQ(cl, library(CLogitTree))
            dev <- parSapply(cl,seeds, one_permutation2,
                          var = variable,
                          exposure = exposure,
                          s = s,
                          kn = knoten,
                          count = count,
                          nvar = nvar,
                          n_levels = n_levels,
                          ordered_values = ordered_values,
                          DM_kov = Z,
                          which_obs = which_obs,
                          splits_evtl = splits_evtl,
                          params = params,
                          dat0 = dat0,
                          mod0 = mod0,
                          n_s = n_s,
                          offset = offset_vec,
                          epsilon = epsilon)
            stopCluster(cl)
          }



          # test decision
          adaption <- sum(!is.na(vars_evtl[[count]][knoten,]))
          crit_val <- quantile(dev,1-(alpha/adaption))
          Tj       <- max(dv[[variable]])
          proof    <- Tj > crit_val
          devs[count]    <- Tj
          crits[count]    <- crit_val
          pvalues[count] <- sum(dev>Tj)/nperm
        } else{
          proof <- TRUE
          devs  <- NULL
          crits <- NULL
          pvalues <- NULL
        }

        if(proof){

          # fit new model
          mod0  <- one_model(variable,exposure,s,knoten,count,split,design_lower,design_upper,params,dat0,offset_vec,epsilon)
          BICs[count+1] <- (-2)*mod0$loglik[2]+log(n)*length(coef(mod0))
          if(count==1 & is.null(exposure)){
            BICs[1] <- (-2)*mod0$loglik[1]+log(n)*length(coef(mod0))
          }
          dat0  <- design[[count+1]] <- data.frame(dat0,design_lower[[variable]][,split,drop=FALSE],design_upper[[variable]][,split,drop=FALSE])
          params_fit[[count+1]] <- names(coefficients(mod0))
# browser()
          # fit final model with penalizedclr
          mm <- model.matrix(mod0, data=dat0)
          if(!is.null(exposure)){
            invisible(capture.output(mod_potential[[count+1]] <- penalized.clr(response=dat0$y, stratum=dat0[,s], penalized=mm[,-c(1,2,ncol(mm))], unpenalized=mm[,2],
                                                                                 offset=offset_vec, lambda=lambda, alpha=10e-20, epsilon = epsilon)))
          } else{
          #  invisible(capture.output(mod_potential[[count+1]] <- penalized.clr(response=dat0$y, stratum=dat0[,s], penalized=mm[,-c(1,2,ncol(mm))], unpenalized=NULL,
          #                                                                       offset=offset_vec, lambda=lambda, alpha=10e-20, epsilon = epsilon)))
            invisible(capture.output(mod_potential[[count+1]] <- penalized.clr(response=dat0$y, stratum=dat0[,s], penalized=mm[,-c(1,ncol(mm))], unpenalized=NULL,
                                                                               lambda=lambda, alpha=10e-20, epsilon = epsilon)))
          }

          # adjust knoten
          if(level>1){
            help_kn4 <- lu(c(),1,level-1,c())
            help_kn5 <- unlist(strsplit(param_old,""))
            help_kn6 <- paste0(help_kn5[which(help_kn5=="_")+1],collapse="")
            knoten2  <- which(help_kn4==help_kn6)
          } else{
            knoten2 <- knoten
          }

          splits[count,"variable"] <- variable
          splits[count,"split"] <- split
          splits[count,"level"] <- level
          splits[count,"node"]  <- knoten2
          splits[count,"threshold"] <- thresholds[[variable]][[split]]
          splits[count,c(6,7,8)] <- c(number,left,right)

          # generiere neue parameter
          params[[count+1]]                             <- params[[count]]
          params[[count+1]]                     <- rep("",length(params[[count]])+1)
          params[[count+1]][c(knoten,knoten+1)] <- param_new
          params[[count+1]][-c(knoten,knoten+1)]<- params[[count]][-knoten]

          # passe splits_evtl an
          n_knots                                                       <- length(params[[count+1]])
          splits_evtl[[count+1]]                                        <- splits_evtl[[count]]
          for(var in 1:nvar){
            splits_evtl[[count+1]][[var]]                       <- matrix(0,nrow=n_knots,ncol=n_s[var])
            splits_evtl[[count+1]][[var]][c(knoten,knoten+1),]  <- matrix(rep(splits_evtl[[count]][[var]][knoten,],2),nrow=2,byrow=T)
            splits_evtl[[count+1]][[var]][-c(knoten,knoten+1),] <- splits_evtl[[count]][[var]][-knoten,]
          }
          splits_evtl[[count+1]][[variable]][knoten,splits_evtl[[count+1]][[variable]][knoten,]>=split] <- NA
          splits_evtl[[count+1]][[variable]][(knoten+1),splits_evtl[[count+1]][[variable]][(knoten+1),]<=split] <- NA

          # passe which_obs an
          which_obs[[count+1]]                               <- which_obs[[count]]
          which_obs[[count+1]]                       <- matrix(0,nrow=n_knots,ncol=n)
          which_obs[[count+1]][c(knoten,knoten+1),]  <- matrix(rep(which_obs[[count]][knoten,],2),nrow=2,byrow=T)
          which_obs[[count+1]][-c(knoten,knoten+1),] <- which_obs[[count]][-knoten,]
          threshh <- ordered_values[[variable]][1:n_s[variable]][split]
          which_obs[[count+1]][knoten,Z[,variable]>threshh] <- NA
          which_obs[[count+1]][(knoten+1),Z[,variable]<=threshh] <- NA

          obs_knoten  <- which(!is.na(which_obs[[count+1]][knoten,]))
          obs_knoten1 <- which(!is.na(which_obs[[count+1]][knoten+1,]))

          # minimal node size constraint
          if(length(obs_knoten)<minnodesize){
            for(var in 1:nvar){splits_evtl[[count+1]][[var]][knoten,] <- rep(NA, n_s[var])}
          }
          if(length(obs_knoten1)<minnodesize){
            for(var in 1:nvar){splits_evtl[[count+1]][[var]][knoten+1,] <- rep(NA, n_s[var])}
          }

          # depth_max
          if(!is.null(depth_max)){
            if(level==depth_max){
              for(var in 1:nvar){splits_evtl[[count+1]][[var]][c(knoten,knoten+1),] <- rep(NA, n_s[var])}
            }
          }

          # table of y
          y_tab[[count+1]]  <- y_tab[[count]]
          y_tab[[count+1]]  <- matrix(0, nrow=n_knots, ncol=2)
          y_tab[[count+1]][knoten,] <- table(factor(y)[obs_knoten])
          y_tab[[count+1]][knoten+1,] <- table(factor(y)[obs_knoten1])
          y_tab[[count+1]][-c(knoten,knoten+1),] <- y_tab[[count]][-knoten,]

          ratio_knoten  <- y_tab[[count+1]][knoten,2]/y_tab[[count+1]][knoten,1]
          ratio_knoten1 <- y_tab[[count+1]][knoten+1,2]/y_tab[[count+1]][knoten+1,1]
          OR_split <- ratio_knoten/ratio_knoten1
          if(!is.na(OR_split)){
            if(OR_split<1){
              splits[count,"association"] <- "-"
            } else{
              splits[count,"association"] <- "+"
            }
          }

          # passe vars_evtl an
          vars_evtl[[count+1]]                     <- vars_evtl[[count]]
          vars_evtl[[count+1]]                     <- matrix(0,nrow=n_knots,ncol=nvar)
          vars_evtl[[count+1]][c(knoten,knoten+1),]  <- matrix(rep(vars_evtl[[count]][knoten,],2),nrow=2,byrow=T)
          vars_evtl[[count+1]][-c(knoten,knoten+1),] <- vars_evtl[[count]][-knoten,]

          vars_evtl[[count+1]][knoten,sapply(1:nvar, function(var) all(is.na(unlist(splits_evtl[[count+1]][[var]][knoten,]))))] <- NA
          vars_evtl[[count+1]][knoten+1,sapply(1:nvar, function(var) all(is.na(unlist(splits_evtl[[count+1]][[var]][knoten+1,]))))] <- NA

          # any split?
          anysplit <- !all(is.na(unlist(splits_evtl[[count+1]])))

          # passe numbers an
          numbers[[count+1]]                                  <- numbers[[count]]
          numbers[[count+1]]                      <- numeric(length=n_knots)
          numbers[[count+1]][c(knoten,knoten+1)]  <- c(left,right)
          numbers[[count+1]][-c(knoten,knoten+1)] <- numbers[[count]][-knoten]

          # trace
          if(print.trace){
            cat(paste0("\n Split"," ",count,";"," ","Variable"," ",variable,"\n"))
          }

          if(length(split)>1){
            warning(paste("Maximum in iteration ",count," not uniquely defined"))
          }

          # erhoehe counter
          count <- count+1

        } else{
          sig <- FALSE
        }
      } else{
        sig <- FALSE
      }
    }
  }
  if(count>1){
    params_fit <- check_names_list(params_fit, params)
  }

  ###################################################################################
  mod_opt         <- mod_potential[[count]]
  params_opt      <- params[[count]]
  params_opt_fit  <- params_fit[[count]]

  if(count==1){
    if(!is.null(exposure)){
      beta_hat <- unlist(mod_opt$penalized)
    } else{
      beta_hat <- NULL
    }
    gamma_hat <- NULL
  }

  if(count>1){
    if(!is.null(exposure)){
      gamma_hat <- c(unlist(mod_opt$penalized),0)
      beta_hat  <- unlist(mod_opt$unpenalized)
      names(gamma_hat) <- params_opt_fit[-1]
      names(beta_hat)  <- params_opt_fit[1]
    } else{
      gamma_hat <- c(unlist(mod_opt$penalized),0)
      beta_hat  <- NULL
      names(gamma_hat) <- params_opt_fit
    }
    gamma_hat <- gamma_hat[params_opt]
  }


  splits$variable  <- names(Z)[splits$variable]


  ## calculate gamma with symmetric side constraints

  if(count>1){
    gamma_hat_sym <- gamma_hat - mean(gamma_hat)
  } else{
    gamma_hat_sym <- NULL
  }

  if(fit){
    to_return <-  list("beta_hat"=beta_hat,
                       "gamma_hat"=gamma_hat,
                       "gamma_hat_sym" = gamma_hat_sym,
                       "splits"=splits,
                       "Z"=Z,
                       "y"=y,
                       "y_tab"=y_tab,
                       "strata" = dat0[,s],
                       "exposure" = dat0[,exposure],
                       "model"=mod_potential,
                       "design"=design,
                       "param"=params,
                       "param_fit"=params_fit,
                       "pvalue"=pvalues,
                       "dev"=devs,
                       "crit"=crits,
                       "BIC"=BICs,
                       "minnodesize" = minnodesize,
                       "minbucket" = minbucket,
                       "depth_max" = depth_max,
                       "lambda" = lambda,
                       "print.trace" = print.trace,
                       "nperm" = nperm,
                       "alpha" = alpha,
                       "perm_test" = perm_test,
                       "mtry" = mtry,
                       "ncores" = ncores,
                       "call"=match.call(),
                       "prunedBIC" = FALSE,
                       "epsilon" = epsilon)
  }else{
    to_return <- list("beta_hat"=beta_hat,
                      "gamma_hat"=gamma_hat,
                      "gamma_hat_sym" = gamma_hat_sym,
                      "splits"=splits,
                      "Z"=Z,
                      "y"=y,
                      "strata" = dat0[,s],
                      "exposure" = dat0[,exposure],
                      "minnodesize" = minnodesize,
                      "minbucket" = minbucket,
                      "depth_max" = depth_max,
                      "lambda" = lambda,
                      "print.trace" = print.trace,
                      "nperm" = nperm,
                      "alpha" = alpha,
                      "perm_test" = perm_test,
                      "mtry" = mtry,
                      "ncores" = ncores,
                      "call"= match.call(),
                      "prunedBIC" = FALSE,
                      "epsilon" = epsilon)
  }

  class(to_return) <- "CLogitTree"
  return(to_return)

}



#' Illustrative data
#'
#' A small simulated dataset for illustration of the CLogitTree method.
#'
#' @name illu.small
#' @docType data
#'
#' @format A data frame with 400 rows and 5 variables
#' \describe{
#'   \item{x}{Exposure variable}
#'   \item{strata}{Strata variable}
#'   \item{Z4}{Explanatory/confounder variable}
#'   \item{Z5}{Explanatory/confounder variable}
#'   \item{y}{Case-control status (1: case; 0: control)}
#' }
NULL


#' Illustrative data
#'
#' A simulated dataset for illustration of the CLogitTree method.
#'
#' @name illu.data
#' @docType data
#'
#' @format A data frame with 1600 rows and 8 variables
#' \describe{
#'   \item{x}{Exposure variable}
#'   \item{strata}{Strata variable}
#'   \item{Z1}{Explanatory/confounder variable}
#'   \item{Z2}{Explanatory/confounder variable}
#'   \item{Z3}{Explanatory/confounder variable}
#'   \item{Z4}{Explanatory/confounder variable}
#'   \item{Z5}{Explanatory/confounder variable}
#'   \item{y}{Case-control status (1: case; 0: control)}
#' }
NULL
