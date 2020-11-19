

CLogitTree <- function(y,
                       X,
                       exposure,
                       s,
                       alpha,
                       nperm,
                       trace=TRUE){
# browser()
  Z    <- X[,!names(X)%in%c(exposure,s)]
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
      Z[,j] <- mod_factors(y,Z[,j])
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
                              "right"=numeric(),stringsAsFactors=FALSE)

  params        <- list()
  vars_evtl     <- list()
  splits_evtl   <- list()
  which_obs     <- list()
  numbers       <- list()
  count         <- 1

  params[[1]]      <- "int"
  which_obs[[1]]   <- matrix(1:n,nrow=1)
  vars_evtl[[1]]   <- nvar
  splits_evtl[[1]] <- lapply(1:nvar, function(var) matrix(1:n_s[var],nrow=1))
  numbers[[1]]     <- 1

  dat0   <- data.frame("y"=y,"int"=rep(1,n),X)
  form0  <- formula(paste0("y~",exposure))
  mod0   <- mod_potential[[1]] <- clogistic(form0, strata = dat0[,s], data = dat0)

  design_lower <- designlists(Z,nvar,n_s,n_levels,ordered_values)[[1]]
  design_upper <- designlists(Z,nvar,n_s,n_levels,ordered_values)[[2]]
  sig      <- TRUE
  anysplit <- TRUE

  while(sig & anysplit){

    # estimate all models
    dv <- lapply(1:nvar,function(var){
              n_knots   <- length(params[[count]])
              deviances <- matrix(rep(0,n_s[var]*n_knots),ncol=n_knots)
              for(kn in 1:n_knots){
                deviances[,kn] <- allmodels(var,exposure,s,kn,count,design_lower,design_upper,splits_evtl,params,dat0,mod0,n_s)
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
    dev <- rep(NA,nperm)

    for(perm in 1:nperm){
      dev[perm] <- one_permutation(variable,exposure,s,knoten,count,nvar,n_levels,ordered_values,
                                   Z,which_obs,splits_evtl,params,dat0,mod0,n_s)
      if(trace){
        cat(".")
      }
    }

    # test decision
    adaption <- vars_evtl[[count]][knoten]
    crit_val <- quantile(dev,1-(alpha/adaption))
    Tj       <- max(dv[[variable]])
    proof    <- Tj > crit_val
    devs[count]    <- Tj
    crits[count]    <- crit_val
    pvalues[count] <- sum(dev>Tj)/nperm

    if(proof){

      # fit new model
      mod0  <- mod_potential[[count+1]] <- one_model(variable,exposure,s,knoten,count,split,design_lower,design_upper,params,dat0)

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

      # any split?
      anysplit <- !all(is.na(unlist(splits_evtl[[count+1]])))

      # passe vars_evtl an
      vars_evtl[[count+1]]                             <- vars_evtl[[count]]
      vars_evtl[[count+1]]                     <- rep(0,n_knots)
      vars_evtl[[count+1]][c(knoten,knoten+1)] <- rep(vars_evtl[[count]][knoten],2)
      vars_evtl[[count+1]][-c(knoten,knoten+1)]<- vars_evtl[[count]][-knoten]

      if(length(which(!is.na(splits_evtl[[count+1]][[variable]][knoten,])))==0){
        vars_evtl[[count+1]][knoten] <- vars_evtl[[count+1]][knoten]-1
      }
      if(length(which(!is.na(splits_evtl[[count+1]][[variable]][knoten+1,])))==0){
        vars_evtl[[count+1]][knoten+1] <- vars_evtl[[count+1]][knoten+1]-1
      }

      # passe which_obs an
      which_obs[[count+1]]                               <- which_obs[[count]]
      which_obs[[count+1]]                       <- matrix(0,nrow=n_knots,ncol=n)
      which_obs[[count+1]][c(knoten,knoten+1),]  <- matrix(rep(which_obs[[count]][knoten,],2),nrow=2,byrow=T)
      which_obs[[count+1]][-c(knoten,knoten+1),] <- which_obs[[count]][-knoten,]
      threshh <- ordered_values[[variable]][1:n_s[variable]][split]
      which_obs[[count+1]][knoten,Z[,variable]>threshh] <- NA
      which_obs[[count+1]][(knoten+1),Z[,variable]<=threshh] <- NA

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

  ###################################################################################
  mod_opt     <- mod_potential[[count]]
  params_opt  <- params[[count]]

  beta_hat <- coefficients(mod_opt)[1]

  if(count>1){
    gamma_hat <- coefficients(mod_opt)[-1]
    gamma_hat[is.na(gamma_hat)] <- 0
    gamma_hat <- gamma_hat[params_opt]
  } else{
    gamma_hat <- NULL
  }

  splits$variable  <- names(Z)[splits$variable]

  to_return <- (list("model"=mod_opt,
                     "beta_hat"=beta_hat,
                     "gamma_hat"=gamma_hat,
                     "splits"=splits,
                     "pvalues"=pvalues,
                     "devs"=devs,
                     "crits"=crits,
                     "Z"=Z,
                     "y"=y,
                     "call"=match.call()))

  class(to_return) <- "CLogitTree"
  return(to_return)

}
