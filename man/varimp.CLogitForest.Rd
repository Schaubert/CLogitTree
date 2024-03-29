% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/varimp.R
\name{varimp.CLogitForest}
\alias{varimp.CLogitForest}
\title{Variable importance for CLogitForest}
\usage{
\method{varimp}{CLogitForest}(model, n.random = 1, oob = TRUE, ...)
}
\arguments{
\item{model}{Original CLogitForest model}

\item{n.random}{Number of randomization iterations used per variable}

\item{oob}{Shall calculation of variable importance values only be based on out-of-bag observations
for each tree}

\item{...}{Further varimp arguments}
}
\value{
Estimated values of variable importance, separately for all variables
}
\description{
Estimates variable importance for objects created by \code{\link{CLogitForest}}.
}
\examples{
\donttest{
data(illu.small)

set.seed(1860)
illu.rf <- CLogitForest(illu.small, response = "y", exposure = "x", s = "strata",
                        ntree = 4, depth_max = 2, tune.mtry = FALSE)

illu.rf

vi <- varimp(illu.rf)
plot(vi)
}
}
\seealso{
\code{\link{CLogitForest}}
}
\author{
Gunther Schauberger: \email{gunther.schauberger@tum.de} \cr
Moritz Berger: \email{moritz.berger@imbie.uni-bonn.de}
}
