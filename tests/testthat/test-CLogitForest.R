test_that("CLogitForest and predict ", {
  n_groups = 50
  samples_per_group = 2
  data <- generate_mixed_dummy_data(n_groups, samples_per_group)  # create a dummy data set with chars, ints and factors
  
  # Test correct output type
  
  model1 <- CLogitForest(data, response = "y", 
                        exposure = "x1", s = "strata", print.trace = FALSE, lambda=1e-20, 
                        ncores = 8, minnodesize = 30, depth_max = 3, minbucket = 10,
                        ntree = 10, mtry = max(floor(sqrt(ncol(data))),1), tune.mtry = FALSE, BIC = FALSE, linear.offset = TRUE)

  model2 <- CLogitForest(data, response = "y", 
                       exposure = "x_factor", s = "strata", print.trace = FALSE, lambda=1e-20, 
                       ncores = 8, minnodesize = 30, depth_max = 3, minbucket = 10,
                       ntree = 10, mtry = max(floor(sqrt(ncol(data))),1), tune.mtry = FALSE, BIC = FALSE, linear.offset = TRUE)

  

  
  models <- list(model1, model2)
  for (model in models){
  # Test attribute structure
  expect_true(is(model, "CLogitForest"))
  expect_true("tree_list" %in% names(model))
  expect_false(anyNA(model$beta_hat))
  expect_equal(length(model$tree_list), 10)  # Should create exactly as many trees as specified
  # Test error handling (needs to be implemented in the first place)
  # expect_error(CLogitForest(data, exposure = NULL, response = "y", s = "strata"), "Invalid exposure, please correct.")
  # expect_error(CLogitForest(data, exposure = "x1", response = NULL, s = "strata"), "Invalid outcome, please correct.")
  
  }
  
  ll1 <- predict(model1, "y", "x1", "strata", newdata = data,type="loglik")
  ll2 <- predict(model2,"y", "x_factor", "strata", newdata = data, type = "loglik")
  lllist <- list(ll1, ll2)
  
  for ( ll in lllist){
    # check if the loglik is always between 0 and 1
    expect_true(all(sapply(ll, function(x) x >= 0)))
    expect_true(all(sapply(ll, function(x) x <= 1)))
    # does it calculcate a loglik for each strata?
    expect_equal(length(ll), n_groups)
    expect_false(anyNA(ll))
  }
  })