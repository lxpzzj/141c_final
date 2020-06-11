#' @import purrr
#' @import Rcpp
#' @import stats
#' @import parallel
#' @import future
#' @importFrom data.table fread
#' @import MASS
#' @importFrom magrittr %>%
#' @aliases blblm.core
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' @title Linear Regression with Little Bag of Bootstraps
#' @description The core function to implement linear Regression with blb.
#' @import RcppArmadillo
#' @param formula The formula to fit.
#' @param data The data from which the regression is fitted
#' @param m The number of bags to divid to
#' @param B The amount of times boostrap will run
#' @export
blblm <- function(formula, data, m = 10, B = 5000) {
  if(typeof(data) == "character"){
    df <- read_list(data)
    data_list <- split_data(df, m)
  }
  else{
    data_list <- split_data(data, m)
  }
  estimates <- map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}
#' @title Linear Regression with Parallel Little Bag of Bootstraps by furrr
#' @description The core function to parallely implement linear Regression with blb by furrr
#' @import future
#' @import furrr
#' @param formula The formula to fit.
#' @param data The data from which the regression is fitted
#' @param m The number of bags to divid to
#' @param B The amount of times boostrap will run
#' @param core The number of core parallel computation will use
#' @export
blblm_furrr <- function(formula, data, m = 10, B = 5000, core = 1) {
  future::plan(future::multiprocess, workers = core)
  if(typeof(data) == "character"){
    df <- read_list(data)
    data_list <- split_data(df, m)
  }
  else{
    data_list <- split_data(data, m)
  }
  estimates <- future_map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}
#' @title Linear Regression with Parallel Little Bag of Bootstraps by parallel package
#' @description The core function to parallely implement linear Regression with blb by parallel
#' @param formula The formula to fit.
#' @param data The data from which the regression is fitted
#' @param m The number of bags to divid to
#' @param B The amount of times boostrap will run
#' @param core The number of core parallel computation will use
#' @export
blblm_parallel <- function(formula, data, m = 10, B = 5000, core = 1) {
  cl <- makeCluster(core)
  clusterExport(cl, c("lm_each_subsample","lm_each_boot","lm1","blbcoef","blbsigma"),envir = environment())
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
                         lm_each_subsample,formula = formula, n = nrow(data), B = B)
  stopCluster(cl)
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' split data into m parts of approximated equal sizes
#' @param data The input data
#' @param m The numbers of subsamples to create
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}

#' read given list of dataset
#' @param data the list of files you want to rwad
read_list <- function(data) {
  pattern = paste0(data, collapse = "|") # Translating name into the regExpress
  # Read all files under working directory at once
  files <- list.files(path = getwd(), pattern = pattern, full.names = T) %>%
    map_df(~fread(.))
  files # Returning the file
}

#' compute the estimates
#' @param formula The linear function to fit
#' @param data The data inputed
#' @param n The numbers of rows to be sampled from
#' @param B The numbers of trails
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}

#' compute the regression estimates for a blb dataset
#' @param formula The linear function to fit
#' @param data The data inputed
#' @param n The numbers of rows to be sampled from
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}

#' estimate the regression estimates based on given the number of repetitions
#' @param formula The linear function to fit
#' @param data The data inputed
#' @param freqs The frequences to boot the data
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}

#' compute the coefficients from fit
#' @param fit The fitted model
blbcoef <- function(fit) {
  coef(fit)
}

#' compute sigma from fit
#' @param fit The fitted model
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' compute sigma from fit
#' @param fit The fitted model
#' @param weights The weights used to fit
#' @param y The response varibale
blbsigma_rcpp <- function(fit, weights,y) {
  p <- length(fit$coefficients)
  y <- y
  e <- fitted(fit) - y
  w <- weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' @title Model Detail
#' @export
#' @param x The fitted model
#' @param ... not use
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

#' @title Get the sigma of the model
#' @export
#' @param object The fitted model
#' @param confidence Whether you want CI or not
#' @param level significant level for computing CI
#' @param ... not use
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}
#' @title Get the coefficiences of the model
#' @export
#' @param object The mode you fitted
#' @param ... not use
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' @title Get the confidence interval of the model
#' @export
#' @param object The model you fitted
#' @param parm The specific parameter you want to use
#' @param level significant level
#' @param ... not use
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(fit$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}
#' @title Predict through the fitted model
#' @export
#' @param object,new_data,confidence,level data
#' @param ... not use
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean_parallel(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}

mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}



