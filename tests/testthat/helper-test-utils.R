# tests/testthat/helper-test-utils.R
generate_dummy_data <- function(n_groups = 50, samples_per_group = 2) {
  set.seed(123)  # For reproducibility
  
  strata <- rep(1:n_groups, each = samples_per_group)
  n <- length(strata)  # Total number of row
  x <- matrix(rnorm(n * 5), n, 5)
  y <- sample(c(0, 1), n, replace = TRUE)  # Binary outcome
  list(x = x, y = y)
}

# A larger dataset with factors
generate_mixed_dummy_data <- function(n_groups = 50, samples_per_group = 2) {
  set.seed(123)  # For reproducibility
  
  strata <- rep(1:n_groups, each = samples_per_group)
  n <- length(strata)  # Total number of row
  
  y <- unlist(lapply(1:n_groups, function(group) rep(c(0, 1), length.out = samples_per_group)))
  
  
  # Numeric data: Continuous values
  x_numeric <- matrix(rnorm(n * 3), n, 3)  # 3 numeric features
  
  # Integer data: Discrete values
  x_integer <- sample(1:10, n, replace = TRUE)  # Random integers between 1 and 10
  
  # Categorical data: Factor
  x_factor <- factor(sample(c("Low", "Low-Medium", "Medium", "Medium-High", "High"), n, replace = TRUE))
  
  
  # Combine all columns into a data frame
  data <- data.frame(
    x1 = x_numeric[, 1],
    x2 = x_numeric[, 2],
    x3 = x_numeric[, 3],
    x_integer = x_integer,
    x_factor = x_factor,
    y = y,
    strata = strata
  )
  
  return(data)
}
