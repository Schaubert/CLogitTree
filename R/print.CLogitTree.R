#' @export
print.CLogitTree <- function(x, ...){

  n.terminal <- length(x$gamma_hat)
  

  cat("Conditional Logistic Regression tree with", n.terminal, "terminal nodes \n \n")

  ### 04.12.2024 JW: Changed the presentation of beta_hat values, added a loop in the else
  if(is.null(x$beta_hat)){
    cat("No separate exposure effect was estimated!","\n","\n")
  } else{
    cat("Exposure effect estimate:\n")
    for (i in seq_along(x$beta_hat)) {
      cat(paste0(names(x$beta_hat)[i], ": "), x$beta_hat[i], "\n", "\n")
    }
  }


  if(!is.null(x$gamma_hat)){

  info <- x$splits
max.level <- max(info$level)

print.frame <- data.frame(threshold = rep(info$threshold[info$level==1],2), variable = rep(info$variable[info$level==1],2),
                          level = rep(1,2), dir = c(1,-1),number = 2:3)


  for(i in 2:max.level){
    splits.i <- info[which(info$level==i),]
    open.nodes <- print.frame[print.frame$level==(i-1),"number"]
    for(j in 1:length(open.nodes)){
      if(open.nodes[j] %in% splits.i$number){
        # right.j <- open.nodes[j] %in% info$right
        j.line <- which(splits.i$number %in% open.nodes[j])
        newline <- rbind(splits.i[j.line,c("threshold", "variable", "level")],
                         splits.i[j.line,c("threshold", "variable", "level")])
        newline <- cbind(newline,c(1,-1),c(splits.i[j,"left"], splits.i[j,"right"]))
        names(newline) <- c("threshold", "variable", "level", "dir", "number")
        i.row <- which(print.frame$number == open.nodes[j])
        print.frame <- insert.row(print.frame, i.row, newline)
      }
    }

  }



    cat("Tree structure:\n")
    for(i in 1:nrow(print.frame)){
      cat(paste0(rep("   ", as.numeric(print.frame[i,"level"])-1)), paste0(i,")"),
          print.frame[i,"variable"], ifelse(print.frame[i,"dir"]>0, "<=", ">"),
          print.frame[i,"threshold"] ,"\n")
    }
  }else{
    cat("Tree without splits!","\n")
  }

}

insert.row <- function(data, i, newline){
  data.lower <- data[1:i,,drop = FALSE]
  if((i+1)<=nrow(data)){
    data.upper <- data[(i+1):nrow(data),,drop = FALSE]
    data.new <- rbind(data.lower, newline, data.upper)
  }else{
    data.new <- rbind(data.lower, newline)
  }

  data.new

}
