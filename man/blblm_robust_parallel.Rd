% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm_robust.R
\name{blblm_robust_parallel}
\alias{blblm_robust_parallel}
\title{Robust Linear Regression with Little Bag of Bootstraps parallely}
\usage{
blblm_robust_parallel(formula, data, m = 10, B = 5000, core = 1)
}
\arguments{
\item{formula}{The formula to fit.}

\item{data}{The data from which the regression is fitted}

\item{m}{The number of bags to divid to}

\item{B}{The amount of times boostrap will run}

\item{core}{The numer of cores to use}
}
\description{
The function to parallely implement robust linear Regression fitting with blb.
}
