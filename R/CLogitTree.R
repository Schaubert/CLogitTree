

CLogitTree <- function(data,
                       response,
                       exposure=NULL,
                       s,
                       alpha,
                       nperm,
                       minnodesize=5,
                       perm_test=TRUE,
                       mtry=NULL,
                       lambda=0,
                       trace=TRUE,
                       fit=TRUE){

  y <- data[, names(data) == response]
  Z    <- data[,!names(data)%in%c(response, exposure,s)]
  n    <- length(y)
  nvar <- ncol(Z)

  if(!is.null(names(Z))){
    var_names <- names(Z)
  } else{
    var_names <- paste0("x",1:nvar)
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
    form0  <- formula(paste0("y~",exposure))
    mod0   <- clogistic(form0, strata = dat0[,s], data = dat0)
    mm <- model.matrix(mod0, data=dat0)
    mod_potential[[1]] <- penalized.clr(response=dat0$y, stratum=dat0[,s], penalized=mm[,2, drop=FALSE], unpenalized=NULL,
                                        lambda=0, alpha=10e-20)
    params_fit[[1]] <- names(coefficients(mod0))
  } else{
    mod0   <- mod_potential[[1]] <- NULL
    params_fit[[1]] <- NULL
  }

  design_lower <- designlists(Z,nvar,n_s,n_levels,ordered_values)[[1]]
  design_upper <- designlists(Z,nvar,n_s,n_levels,ordered_values)[[2]]
  sig      <- TRUE
  anysplit <- TRUE

  while(sig & anysplit){

    n_knots <- length(params[[count]])
    if(!is.null(mtry)){
      tryvars <- sapply(1:n_knots, function(kn) sample(1:nvar, mtry))
    } else{
      tryvars <- sapply(1:n_knots, function(kn) 1:nvar)
    }

    # estimate all models
    dv <- lapply(1:nvar,function(var){
              deviances <- matrix(rep(0,n_s[var]*n_knots),ncol=n_knots)
              which_knots <- which(!is.na(vars_evtl[[count]][,var]))
              for(kn in which_knots){
                if(var %in% tryvars[,kn]){
                  deviances[,kn] <- allmodels(var,exposure,s,kn,count,design_lower,design_upper,splits_evtl,params,dat0,mod0,n_s)
                }
              }
              return(deviances)
          })


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
      dev <- rep(NA,nperm)

      for(perm in 1:nperm){
        dev[perm] <- one_permutation(variable,exposure,s,knoten,count,nvar,n_levels,ordered_values,
                                     Z,which_obs,splits_evtl,params,dat0,mod0,n_s)
        if(trace){
          cat(".")
        }
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
      mod0  <- one_model(variable,exposure,s,knoten,count,split,design_lower,design_upper,params,dat0)
      dat0  <- design[[count+1]] <- data.frame(dat0,design_lower[[variable]][,split,drop=FALSE],design_upper[[variable]][,split,drop=FALSE])
      params_fit[[count+1]] <- names(coefficients(mod0))

      # fit final model with penalizedclr
      mm <- model.matrix(mod0, data=dat0)
      if(!is.null(exposure)){
        mod_potential[[count+1]] <- penalized.clr(response=dat0$y, stratum=dat0[,s], penalized=mm[,-c(1,2,ncol(mm))], unpenalized=mm[,2],
                                                  lambda=lambda, alpha=10e-20)
      } else{
        mod_potential[[count+1]] <- penalized.clr(response=dat0$y, stratum=dat0[,s], penalized=mm[,-c(1,ncol(mm))], unpenalized=NULL,
                                                  lambda=lambda, alpha=10e-20)
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

      # table of y
      y_tab[[count+1]]  <- y_tab[[count]]
      y_tab[[count+1]]  <- matrix(0, nrow=n_knots, ncol=2)
      y_tab[[count+1]][knoten,] <- table(factor(y)[obs_knoten])
      y_tab[[count+1]][knoten+1,] <- table(factor(y)[obs_knoten1])
      y_tab[[count+1]][-c(knoten,knoten+1),] <- y_tab[[count]][-knoten,]

      ratio_knoten  <- y_tab[[count+1]][knoten,2]/y_tab[[count+1]][knoten,1]
      ratio_knoten1 <- y_tab[[count+1]][knoten+1,2]/y_tab[[count+1]][knoten+1,1]
      OR_split <- ratio_knoten/ratio_knoten1
      if(OR_split<1){
        splits[count,"association"] <- "-"
      } else{
        splits[count,"association"] <- "+"
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
      if(trace){
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
  }
  params_fit <- check_names_list(params_fit, params)

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

    gamma_hat_sym <- gamma_hat - mean(gamma_hat)


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
                       "minnodesize" = minnodesize,
                       "lambda" = lambda,
                       "trace" = trace,
                       "nperm" = nperm,
                       "alpha" = alpha,
                       "perm_test" = perm_test,
                       "mtry" = mtry,
                       "call"=match.call())
  }else{
    to_return <- list("beta_hat"=beta_hat,
                      "gamma_hat"=gamma_hat,
                      "gamma_hat_sym" = gamma_hat_sym,
                      "splits"=splits,
                      "Z"=Z,
                      "y"=y,
                      "y_tab"=y_tab,
                      "minnodesize" = minnodesize,
                      "lambda" = lambda,
                      "trace" = trace,
                      "nperm" = nperm,
                      "alpha" = alpha,
                      "perm_test" = perm_test,
                      "mtry" = mtry,
                      "call"= match.call())
  }

  class(to_return) <- "CLogitTree"
  return(to_return)

}
