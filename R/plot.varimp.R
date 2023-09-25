#' Plot function for variable importance of \code{\link{CLogitForest}}
#'
#' Plots barplots which represent the variable importance in \code{\link{CLogitForest}} objects calculated via \code{\link{varimp.CLogitForest}}.
#'
#' @param x object containing variable importance results for a \code{\link{CLogitForest}} object
#' @param sort Shall the variables be plotted sorted according to their variable importance values?
#' @param horiz See \code{\link{barplot}}
#' @param las See \code{\link{barplot}}
#' @param col See \code{\link{barplot}}
#' @param border See \code{\link{barplot}}
#' @param ... further \code{\link{barplot}} arguments
#' @return No return, used for side effects
#' @author Gunther Schauberger: \email{gunther.schauberger@@tum.de} \cr
#' Moritz Berger: \email{moritz.berger@@imbie.uni-bonn.de}
#' @seealso \code{\link{CLogitForest}}, \code{\link{varimp.CLogitForest}}
#' @examples
#' data(illu.small)
#'
#' @export
plot.varimp.CLF <- function(x,
                            sort = TRUE,
                            horiz = TRUE,
                            las = 1,
                            col = "deepskyblue4",
                            border = FALSE,
                            ...){

  if(sort){x <- sort(x)}

  barplot(c(x),  horiz = horiz, las = las, col = col, border = border, ...)


}
