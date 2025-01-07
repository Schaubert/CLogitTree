#' @export
print.CLogitForest <- function(x, ...){

  cat("CLogitForest consisting of", x$ntree, "trees \n \n ")

  ### JW 05.12.2025 Adjusted the presentation of beta_hat
  if(is.null(x$beta_hat)){
    cat("No separate exposure effect was estimated!","\n","\n")
  }
  else{
    cat("Average exposure effect of exposure variable \n")
    for (i in seq_along(x$beta_hat)) {
      cat(paste0(names(x$beta_hat)[i], ": "), x$beta_hat[i], "\n", "\n")
    }
    cat("Median exposure effect of exposure variable \n")
    for (i in seq_along(x$beta_median_hat)) {
      cat(paste0(names(x$beta_median_hat)[i], ": "), x$beta_median_hat[i], "\n", "\n")
    }
  }
}
