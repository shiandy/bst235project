% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{adaptive_lasso}
\alias{adaptive_lasso}
\title{Fit adaptive lasso model}
\usage{
adaptive_lasso(xs, ys, gam = 1, nfolds = 5)
}
\arguments{
\item{xs}{Matrix of predictors. Should be standardized and centered
to have mean 0, variance 1}

\item{ys}{Matrix or vector of observations. Does not need to be
centered.}

\item{gam}{Power for penalty weights.}

\item{nfolds}{Number of cross-validation folds.}
}
\value{
A cv.glmnet object
}
\description{
Fit adaptive lasso model with intercept. The penalty for each beta_j
coefficient in the model is weighted by 1 / |beta_ini_j|^gam, where
beta_ini_j is an initial estimate of beta, here obtained through OLS.
Then, a LASSO is run. The optimal lambda is selected through
cross-validation.
}
