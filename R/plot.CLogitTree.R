#' Plot function for CLogitTree
#'
#' Plots trees estimated by \code{\link{CLogitTree}}.
#'
#' @param x CLogitTree object
#' @param symmetric Shall parameter estimates with symmetric side constraints be used?
#' @param precision Number of decimals of parameter estimates
#' @param ellipse_a Controls width of ellipse containing node information
#' @param ellipse_b Controls height of ellipse containing node information
#' @param ellipse_x Controls location on x-axis of ellipse containing node information
#' @param ellipse_y Controls location on y-axis of ellipse containing node information
#' @param branch_adj Vertical adjustment of branch labels
#' @param cex.lines Size of lines
#' @param cex.branches Size of branch labels
#' @param cex.coefs Size of parameter estimates
#' @param cex.main Size of title
#' @param cex.numbers Size of numbers of observations
#' @param info_inner Shall +/- signs be added to split points?
#' @param info_n Shall numbers of observations be displayed?
#' @param draw_numbers Shall internally used node numbers be displayed?
#' @param title Main title
#' @param ... further plot arguments
#' @return No return, used for side effects
#' @author Gunther Schauberger: \email{gunther.schauberger@@tum.de} \cr
#' Moritz Berger: \email{moritz.berger@@imbie.uni-bonn.de}
#' @seealso \code{\link{CLogitTree}}, \code{\link{bootci.CLogitTree}}, \code{\link{prune.CLogitTree}}
#' @examples
#' data(illu.small)
#'
#' set.seed(1860)
#' illu.tree <- CLogitTree(illu.small, response = "y", exposure = "x", s = "strata",
#'                         alpha = 0.05, nperm = 20, print.trace = FALSE)
#'
#' plot(illu.tree)
#'
#' @export
plot.CLogitTree <- function(x,
                            symmetric = TRUE,
                            precision = 3,
                            ellipse_a=0.8,
                            ellipse_b=0.2,
                            ellipse_x=0,
                            ellipse_y=0.2,
                            branch_adj=0,
                            cex.lines=2,
                            cex.branches=1,
                            cex.coefs=1,
                            cex.main=1,
                            cex.numbers=1,
                            info_inner=TRUE,
                            info_n=TRUE,
                            draw_numbers=FALSE,
                            title=NULL,
                            ...){

  if(is.null(x$splits)){
    cat("There is no plot available")
  } else{

    X        <- x$Z
    info     <- x$splits

    if(nrow(info)==0){
      cat("There is no tree to plot.\n")
    } else{


      if(symmetric){
        coefs_hat <- x$gamma_hat_sym
      }else{
        coefs_hat <- x$gamma_hat
      }

      endnodes      <- list()
      endnodes[[1]] <- 1
      for(j in 1:nrow(info)){
        endnodes[[j+1]] <- numeric(length=(j+1))
        what <- c(info[j,"left"],info[j,"right"])
        delete <- info[j,"number"]
        where  <- which(endnodes[[j]]==delete)
        endnodes[[j+1]][c(where,where+1)] <- what
        endnodes[[j+1]][-c(where,where+1)] <- endnodes[[j]][-where]
      }
      endnodes <- endnodes[[length(endnodes)]]
      dir      <- sapply(1:length(endnodes), function(j){ifelse(endnodes[j] %in% info[,7], "l","r")})

      n_levels <- length(unique(info[,"level"]))

      hilfspunkte <- list()
      hilfspunkte[[1]] <- matrix(NA,nrow=2^n_levels,ncol=2)
      hilfspunkte[[1]][,1] <- 2^n_levels
      hilfspunkte[[1]][,2] <- rep(n_levels+1,2^n_levels)

      steps <- 2^((n_levels:1-1))

      for(i in 1:n_levels){

        hilfspunkte[[i+1]] <- hilfspunkte[[i]]
        hilfspunkte[[i+1]][,2] <- rep(n_levels+1-i,2^n_levels)

        help  <- c(-steps[i],steps[i])
        help1 <- rep(help,each=steps[i])
        help2 <- rep(help1,length=2^n_levels)
        hilfspunkte[[i+1]][,1] <- hilfspunkte[[i]][,1]+help2

        which_knots <- info[info[,"level"]==i,"node"]
        help3 <- seq(1,2^n_levels)
        help4 <- split(help3,rep(1:2^(i-1),each=2^n_levels/2^(i-1)))
        help5 <- unlist(lapply(which_knots, function(j) help4[[j]]))
        hilfspunkte[[i+1]][-help5,] <- hilfspunkte[[i]][-help5,]

      }


      plot.new(...)
      plot.window(ylim=c(0.5,n_levels+1),xlim=c(0,2^(n_levels+1)))
      # rect(0,0.5,2^(n_levels+1),n_levels+1, border = grey(0.9),col = grey(0.9))


      for(j in length(hilfspunkte):2){
        for(i in 1:(2^n_levels)){
          lines(c(hilfspunkte[[j-1]][i,1],hilfspunkte[[j]][i,1]),c(hilfspunkte[[j-1]][i,2],hilfspunkte[[j]][i,2]),
                lwd=cex.lines)
        }
      }
      if(is.null(title)){
        title <- "f(Z)"
      }
      title(title,cex.main=cex.main)

      # Fuege Schaetzer in den Knoten hinzu
      betas_hat <- format(round(coefs_hat, precision),nsmall=precision)
      y_tab     <- x$y_tab[[length(x$y_tab)]]
      points_betas <- unique(hilfspunkte[[n_levels+1]])
      for(i in 1:length(betas_hat)){
        if(dir[i]=="l"){
          fac <- -1
        } else{
          fac <- 1
        }
        draw.ellipse(x=points_betas[i,1]+fac*ellipse_x,y=points_betas[i,2]-ellipse_y,a=ellipse_a,b=ellipse_b,lwd=cex.lines,col=grey(0.8))
        if(info_n){
          text(points_betas[i,1]+fac*ellipse_x,points_betas[i,2]-ellipse_y,
               bquote(atop(.(betas_hat[i]), n[0]*"="*.(y_tab[i,1])~~n[1]*"="*.(y_tab[i,2]))),cex=cex.coefs)
        } else{
          text(points_betas[i,1]+fac*ellipse_x,points_betas[i,2]-ellipse_y,betas_hat[i],cex=cex.coefs)
        }
      }

      # Fuege Knotennummern hinzu
      if(draw_numbers){
        for(i in 1:length(betas_hat)){
          if(dir[i]=="l"){
            fac <- -1
          } else{
            fac <- 1
          }
          rect(points_betas[i,1]+fac*ellipse_x-max(0.2,ellipse_x/3),points_betas[i,2]+ellipse_b-max(0.05,ellipse_b/3),points_betas[i,1]+fac*ellipse_x+max(0.2,ellipse_x/3),points_betas[i,2]+ellipse_b+max(0.05,ellipse_b/3),col=grey(0.9),lwd=cex.numbers)
          text(points_betas[i,1]+fac*ellipse_x,points_betas[i,2]+ellipse_b,endnodes[i],cex=cex.numbers)
        }
      }

      # Fuege weitere Beschriftung hinzu
      for(i in 1:nrow(info)){
        help4 <- split(help3,rep(1:2^(info[i,"level"]-1),each=2^n_levels/2^(info[i,"level"]-1)))[[info[i,"node"]]]
        point_var <- unique(hilfspunkte[[info[i,"level"]]][help4,])
        points(point_var[1],point_var[2],cex=cex.lines-1,pch=19)
        if(info_inner){
          if(!is.na(info[i,"association"])){
            if(info[i,"association"]=="-"){
              text(point_var[1]-0.5, point_var[2], "-", cex=cex.lines)
              text(point_var[1]+0.5, point_var[2], "+", cex=cex.lines)
            } else{
              text(point_var[1]-0.5, point_var[2], "+", cex=cex.lines)
              text(point_var[1]+0.5, point_var[2], "-", cex=cex.lines)
            }
          }
        }
        point_left  <- c(point_var[1]-steps[info[i,"level"]]-branch_adj,point_var[2]-0.5)
        point_right <- c(point_var[1]+steps[info[i,"level"]]+branch_adj,point_var[2]-0.5)
        var   <- info[i,"variable"]
        thres <- info[i,"threshold"]
        sort_values <- unique(sort(X[,var]))
        if(thres==min(sort_values)){
          text(point_left[1],point_left[2],paste0(var,"=",round(thres,2)),cex=cex.branches,adj=c(1,0))
        } else{
          text(point_left[1],point_left[2],paste0(var,"<=",round(thres,2)),cex=cex.branches,adj=c(1,0))
        }
        if(thres==max(sort_values[-length(sort_values)])){
          text(point_right[1],point_right[2],paste0(var,"=",round(max(sort_values),2)),cex=cex.branches,adj=c(0,0))
        } else{
          text(point_right[1],point_right[2],paste0(var,">",round(thres,2)),cex=cex.branches,adj=c(0,0))
        }
      }
    }
  }
}
