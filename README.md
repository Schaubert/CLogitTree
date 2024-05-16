# CLogitTree

CLogitTree is a package which provides machine-learning methods for the analysis of matched case-control studies. - CLogitTree() is a function to estimate decision trees for matched case-control studies - CLogitForest() is a function to estimate random forests for matched case-control studies

For further information on CLogitTree() see <https://doi.org/10.1002/sim.9637>

# Installation

The package can be installed in R via

``` r
install.packages("devtools")
devtools::install_github("Schaubert/ClogitTree")
```

# Example CLogitTree

``` r
library(CLogitTree)
data(illu.small)

# Example #1: Split selection via permutation trees
set.seed(1860)
illu.tree <- CLogitTree(illu.small, response = "y", 
                        exposure = "x", s = "strata", 
                        alpha = 0.05, nperm = 20)

plot(illu.tree)

# Example #2: Pruning via BIC
 set.seed(1860)
 illu.tree <- CLogitTree(illu.small, response = "y", 
                         exposure = "x", s = "strata", 
                         perm_test = FALSE, depth_max=4)
 
 illu.tree <- pruneBIC(illu.tree)

plot(illu.tree)
```

# Example CLogitForest

``` r
library(CLogitTree)
data(illu.small)

set.seed(1860)
illu.rf <- CLogitForest(illu.small, response = "y", 
                        exposure = "x", s = "strata", 
                        ntree = 4, depth_max = 2, 
                        tune.mtry = FALSE)

illu.rf

vi <- varimp(illu.rf)
plot(vi)
```
