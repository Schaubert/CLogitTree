#' @export
print.CLogitForest <- function(x, ...){

  cat("CLogitForest consisting of", x$ntree, "trees \n \n ")

  cat("Average exposure effect of exposure variable", x$exposure, "\n")
  print(x$beta_hat, ...)


  cat("Median exposure effect of exposure variable", x$exposure, "\n")
  print(x$beta_median_hat, ...)
}
