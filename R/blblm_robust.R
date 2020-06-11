#' @title Robust Linear Regression with Little Bag of Bootstraps
#' @description The function to implement robust linear Regression fitting with blb.
#' @import purrr
#' @param formula The formula to fit.
#' @param data The data from which the regression is fitted
#' @param m The number of bags to divid to
#' @param B The amount of times boostrap will run
#' @export
blblm_robust <- function(formula, data, m = 10, B = 5000) {
  if(typeof(data) == "character"){
    df <- read_list(data)
    data_list <- split_data(df, m)
  }
  else{
    data_list <- split_data(data, m)
  }
  estimates <- map(
    data_list,
    ~ lm_robust_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

lm_robust_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}

lm_robust_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}

lm1_robust <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- rlm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' @title Robust Linear Regression with Little Bag of Bootstraps parallely
#' @description The function to parallely implement robust linear Regression fitting with blb.
#' @import purrr
#' @param formula The formula to fit.
#' @param data The data from which the regression is fitted
#' @param m The number of bags to divid to
#' @param B The amount of times boostrap will run
#' @param core The numer of cores to use
#' @export
blblm_robust_parallel <- function(formula, data, m = 10, B = 5000, core = 1) {
  cl <- makeCluster(core)
  clusterExport(cl, c("lm_robust_each_boot","lm1_robust","blbcoef","blbsigma"),envir = environment())
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
                           lm_robust_each_subsample,formula = formula, n = nrow(data), B = B)
  stopCluster(cl)
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


