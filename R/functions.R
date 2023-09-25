# adjust node
lu <- function(last=last, cd=cd, d=d, erg){

  s <- c("l","u")

  for(i in 1:length(s)){
    s_new <- paste0(last,s[i])
    if(cd==d){
      erg[length(erg)+1] <- s_new
    } else{
      erg <- lu(s_new,cd+1,d,erg)
    }
  }
  return(erg)
}

# modify factors
mod_factors <- function(y, x){
  tab <- table(x,y)
  tab <- tab[rowSums(tab)>0,]
  nx  <- rowSums(tab)
  ptab <- tab/nx
  pp  <- colMeans(ptab, na.rm=T)
  ptabc <- t(apply(ptab,1, function(x) x-pp))
  sig   <- matrix(0, nrow=ncol(tab), ncol=ncol(tab))
  for(j in 1:nrow(tab)){
    sig <- sig+nx[j]*ptabc[j,]%*%t(ptabc[j,])
  }
  sig <- sig/(length(y)-1)
  v <- eigen(sig)$vectors[,1]
  sa <- apply(ptabc, 1, function(x) v%*%x)
  xnew <- as.numeric(x)
  xnew <- sapply(1:length(xnew), function(j) which(order(sa)==xnew[j]))
  return(xnew)
}

# compute ordered values
ord_values <- function(x){
  if(!all((x - round(x)) == 0) || length(unique(x))>50){
    ret <- quantile(x,seq(0.05,1,by=0.05))
  } else{
    ret <- unique(sort(x))
  }
  return(ret)
}

thresh <- function(ordered_values){
  ret <- ordered_values[-length(ordered_values)]
  return(ret)
}


# functions to build design
design_one  <- function(x,threshold,upper){
  if(upper){
    ret <- ifelse(x > threshold,1,0)
  } else{
    ret <- ifelse(x > threshold,0,1)
  }
  return(ret)
}

design <- function(x,thresholds,upper){
  ret <- sapply(thresholds, function(j) design_one(x,j,upper))
  return(ret)
}

designlist <- function(X,vars,label,thresholds,var_names,upper=TRUE){
  ret <- lapply(vars, function(j) {
    ret <- design(X[,j],thresholds[[j]],upper)
    colnames(ret) <- paste0("s",which(var_names==j),1:length(thresholds[[j]]),ifelse(upper,"_u","_l"),label)
    return(ret)})
  names(ret) <- vars
  return(ret)
}

designlists <- function(DM_kov,nvar,n_s,n_levels,ordered_values){

  # generate design matrices for tree
  v <- lapply(1:nvar,function(j) 1:(n_levels[j]-1))
  w <- lapply(1:nvar, function(j) rep(paste0("s",j),n_s[j]))

  design_upper <- lapply(1:nvar, function(j){
    design_matrix <- sapply(1:(n_levels[j]-1),function(k) { ifelse(DM_kov[,j] > ordered_values[[j]][k],1,0)})
    colnames(design_matrix) <- paste0(w[[j]],v[[j]],"_u")
    design_matrix
  })
  design_lower <- lapply(1:nvar, function(j){
    design_matrix <- abs(design_upper[[j]]-1)
    colnames(design_matrix) <- paste0(w[[j]],v[[j]],"_l")
    design_matrix
  })
  return(list(design_lower,design_upper))
}

one_model <- function(var,exposure,s,kn,count,j,design_lower,design_upper,params,dat0,offset_vec, epsilon){

  dat   <- data.frame(dat0,design_lower[[var]][,j,drop=FALSE],design_upper[[var]][,j,drop=FALSE])

  help1 <- params[[count]]
  help2 <- help1[-kn]
  help3 <- paste(help1[kn],c(colnames(design_lower[[var]])[j],colnames(design_upper[[var]])[j]),sep=":")
  help4 <- paste(c(exposure,help2,help3), collapse="+")
  help5 <- formula(paste0("y~", help4, "+ offset(offset_vec)"))
  mod   <- clogistic(help5, strata=dat[,s], data=dat, eps = epsilon, toler.chol = epsilon/10)
  return(mod)
}

allmodels <- function(var,exposure,s,kn,count,design_lower,design_upper,splits_evtl,params,dat0,mod0,n_s,offset_vec,epsilon){
# browser()
  deviances <- rep(0,n_s[var])
  splits_aktuell <- splits_evtl[[count]][[var]][kn,]
  splits_aktuell <- splits_aktuell[!is.na(splits_aktuell)]

  if(length(splits_aktuell)>0){
    for(j in splits_aktuell){
      mod <- try(one_model(var,exposure,s,kn,count,j,design_lower,design_upper,params,dat0,offset_vec,epsilon), silent=T)
      if(!inherits(mod, "try-error")){
        if(is.null(mod0)){
          LRstat <- (-2)*mod$loglik[1] - (-2)*mod$loglik[2]
        } else{
          LRstat <- (-2)*mod0$loglik[2] - (-2)*mod$loglik[2]
        }
        if(!is.na(LRstat)){
          deviances[j] <- LRstat
        }
      }
    }
  }
  return(deviances)
}

one_permutation <- function(var,exposure,s,kn,count,nvar,n_levels,ordered_values,
                            DM_kov,which_obs,splits_evtl,params,dat0,mod0,n_s,offset,epsilon){

  obs_aktuell <- which_obs[[count]][kn,]
  obs_aktuell <- obs_aktuell[!is.na(obs_aktuell)]
  DM_kov_perm <- DM_kov
  DM_kov_perm[obs_aktuell,var] <- sample(DM_kov_perm[obs_aktuell,var],length(obs_aktuell))

  design_upper_perm      <- designlists(DM_kov_perm,nvar,n_s,n_levels,ordered_values)[[1]]
  design_lower_perm      <- designlists(DM_kov_perm,nvar,n_s,n_levels,ordered_values)[[2]]

  dv_perm <- allmodels(var,exposure,s,kn,count,design_lower_perm,design_upper_perm,splits_evtl,params,dat0,mod0,n_s,offset,epsilon)

  return(max(dv_perm))

}

one_permutation2 <- function(seed, var, exposure,s,kn,count,nvar,n_levels,ordered_values,
                            DM_kov,which_obs,splits_evtl,params,dat0,mod0,n_s,offset,epsilon){
# browser()
  set.seed(seed)
  obs_aktuell <- which_obs[[count]][kn,]
  obs_aktuell <- obs_aktuell[!is.na(obs_aktuell)]
  DM_kov_perm <- DM_kov
  DM_kov_perm[obs_aktuell,var] <- sample(DM_kov_perm[obs_aktuell,var],length(obs_aktuell))

  design_upper_perm      <- designlists(DM_kov_perm,nvar,n_s,n_levels,ordered_values)[[1]]
  design_lower_perm      <- designlists(DM_kov_perm,nvar,n_s,n_levels,ordered_values)[[2]]

  dv_perm <- allmodels(var,exposure,s,kn,count,design_lower_perm,design_upper_perm,splits_evtl,params,dat0,mod0,n_s,offset,epsilon)

  return(max(dv_perm))

}

# functions to fix names
check_names <- function(model_names, param_names){

  params_c  <- unlist(param_names)
  names_mod <- unlist(model_names)
  if(!all(params_c %in% names_mod)){
    params_split <- strsplit(params_c,":")
    names_split  <- strsplit(names_mod,":")
    for(k in 1:length(params_split)){
      found <- FALSE
      kk   <- 0
      while(!found){
        kk <- kk+1
        found <- all(params_split[[k]] %in% names_split[[kk]])
      }
      model_names[kk] <- params_c[k]
    }
  }
  return(model_names)
}

check_names_list <- function(model_names_list, param_names_list){
  for(j in 2:length(model_names_list)){
    model_names_list[[j]] <- check_names(model_names_list[[j]],param_names_list[[j]])
  }
  return(model_names_list)
}


one_boot_fun <- function(seed, index, Xboot, alpha, nperm, minnodesize, minbucket, depth_max,
                         perm_test, mtry, lambda, print.trace, fit, prunedBIC, epsilon){
  set.seed(seed)

  index2 <- sample(index, size = length(index), replace = TRUE)

  X2 <- c()
  for(i in 1:length(index2)){
    Xnew <- Xboot[Xboot$Xstrata == index2[i],]
    Xnew$Xstrata <- rep(i, nrow(Xnew))
    X2 <- rbind(X2, Xnew)
  }

  ret <- CLogitTree(data = X2,
                    response = "Xy",
                    exposure="Xexpo",
                    s = "Xstrata",
                    alpha = alpha,
                    nperm = nperm,
                    minnodesize=minnodesize,
                    minbucket = minbucket,
                    depth_max = depth_max,
                    perm_test=perm_test,
                    mtry=mtry,
                    lambda=lambda,
                    print.trace=print.trace,
                    fit=TRUE,
                    ncores = 1,
                    epsilon = epsilon)

  if(prunedBIC){
    ret <- pruneBIC(ret)
  }

  return(ret$beta_hat)

}



one_boot_forest <- function(seed, index, Xboot, minnodesize, minbucket, depth_max,
                         ntree, mtry, subsample.632, lambda, print.trace, fit, response,
                         exposure, s, BIC, linear.offset,epsilon){
  set.seed(seed)


  index2 <- sample(index, size = length(index), replace = TRUE)

  X2 <- c()
  for(i in 1:length(index2)){
    Xnew <- Xboot[Xboot[,names(Xboot) == s] == index2[i],]
    Xnew[,names(Xnew) == s] <- rep(i, nrow(Xnew))
    X2 <- rbind(X2, Xnew)
  }

cat("Start new forest on bootstrap:\n")

  ret <- CLogitForest(data = X2,
                    response = response,
                    exposure=exposure,
                    s = s,
                    ntree = ntree,
                    mtry = mtry,
                    subsample.632 = subsample.632,
                    minnodesize=minnodesize,
                    minbucket = minbucket,
                    depth_max = depth_max,
                    lambda = lambda,
                    print.trace = print.trace,
                    fit = TRUE,
                    ncores = 1,
                    BIC = BIC,
                    # offset = offset,
                    linear.offset = linear.offset,
                    tune.mtry=FALSE,
                    epsilon = epsilon)

  return(ret$beta_hat)

}
