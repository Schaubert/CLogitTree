test_that("CLogitTree and predict", {
  n_groups = 50
  samples_per_group = 2
  data <- generate_mixed_dummy_data(n_groups, samples_per_group)  # create a dummy data set with chars, ints and factors
  
  # Test correct output type
  
  model <- CLogitTree(data, response = "y", 
                                 exposure = "x1", s = "strata", alpha = 0.1, 
                                 perm_test = FALSE, print.trace = FALSE, lambda=1e-30, 
                                 ncores = 4, minnodesize = 15, depth_max = 4, nperm = 0, minbucket = 5)
  modelBIC <- pruneBIC(model)
  
  model1 <- CLogitTree(data, response = "y", 
                      exposure = "x_factor", s = "strata", alpha = 0.1, 
                      perm_test = FALSE, print.trace = FALSE, lambda=1e-30, 
                      ncores = 4, minnodesize = 15, depth_max = 4, nperm = 0, minbucket = 5)
  model1BIC <- pruneBIC(model1)
  
  
  models <- list(modelBIC, model1BIC)
  
  lltreeBIC <- predict(model1BIC, "y", "x_factor", "strata", newdata = data, type="loglik")
  lltree1BIC <- predict(modelBIC,"y", "x1", "strata", newdata = data, type = "loglik")
  lllist <- list(lltreeBIC, lltree1BIC)
  
  for ( ll in lllist){
  # check if the loglik is always between 0 and 1
  expect_true(all(sapply(ll, function(x) x >= 0)))
  expect_true(all(sapply(ll, function(x) x <= 1)))
  # does it calculcate a loglik for each strata?
  expect_equal(length(ll), n_groups)
  expect_false(anyNA(ll))
  }
  
  for (model in models){
  expect_true(is(model, "CLogitTree")) # test if a CLogitTree object is created

  # Test if an effect was calculated
  expect_false(anyNA(model$beta_hat))
  }

  # Test error handling (needs to be implemented in the first place)
  # expect_error(CLogitTree(data, exposure = NULL, response = "y", s = "strata"), "Invalid exposure, please correct.")
  # expect_error(CLogitTree(data, exposure = "x1", response = NULL, s = "strata"), "Invalid outcome, please correct.")
  })
