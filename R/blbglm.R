#' @title Generalized Linear Regression with Little Bag of Bootstraps
#' @description This function implement the GLM model by blb method
#' @import purrr
#' @import Rcpp
#' @import stats
#' @importFrom magrittr %>%
#' @param formula an object of class "formula" (or one that can be coerced to that class).
#' @param family a description of the error distribution and link function to be used in the model.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param m The amounts of subsamples to create
#' @param B The times of boostrapy
#' #' @export
blbglm <- function(formula, data, family = 'gaussian', m = 10, B = 5000) {
  if(typeof(data) == "character"){
    df <- read_list(data)
    data_list <- split_data(df, m)
  }
  else{
    data_list <- split_data(data, m)
  }
  estimates <- map(
    data_list,
    ~ glm_each_subsample(formula = formula, family = family, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}
#' @title Generalized Linear Regression with Little Bag of Bootstraps parallely
#' @description This function implement the GLM model by blb method with parallel computation
#' @import purrr
#' @import Rcpp
#' @import stats
#' @import parallel
#' @importFrom magrittr %>%
#' @param formula an object of class "formula" (or one that can be coerced to that class).
#' @param family a description of the error distribution and link function to be used in the model.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.
#' @param m The amounts of subsamples to create
#' @param B The times of boostrapy
#' @param core The numbers of cores to use for parallel computation
#' @export
blbglm_prallel <- function(formula, data, family = 'gaussian', m = 10, B = 5000, core = 1) {
  cl <- makeCluster(core)
  clusterExport(cl, c("glm_each_boot","glm_base","blbcoef","blbsigma"),envir = environment())
  if(typeof(data) == "character"){
    df <- read_list(data)
    data_list <- split_data(df, m)
  }
  else{
    data_list <- split_data(data, m)
  }
  estimates <- parLapplyLB(cl,
                           data_list,
                           chunk.size = 2,
                           glm_each_subsample,formula = formula, family = family, n = nrow(data), B = B)
  stopCluster(cl)
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

glm_each_subsample <- function(formula, family, data, n, B) {
  replicate(B, glm_each_boot(formula, family = family, data, n), simplify = FALSE)
}

glm_each_boot <- function(formula, family, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm_base(formula, data, freqs, family = family)
}

glm_base <- function(formula, family, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- glm(formula, family = family, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}



