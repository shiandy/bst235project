#' Simulate violations of link function
#'
#' Simulate violations of link function. The true model is
#' E(Y_i | X_i) = alpha_0 + Phi(beta_0^T X_i),
#' where Phi is the Normal CDF, but we fit the model
#' E(Y_i | X_i) = alpha + beta^T X
#' using adaptive LASSO.
#'
#' We calculate the errors using a large validation dataset. The errors
#' considered are
#' \enumerate{
#' \item E((Y_i - hat(alpha) - hat(beta)^T X_i)^2), where hat(alpha) and
#' hat(beta) are from adaptive LASSO.
#' \item E((Y_i - tilde(m)(X_i)^2), where tilde(m)(X_i) corresponds to
#' the conditional mean from fitting the true mode.
#' \item E((Y_i - hat(m)(X_i))^2), where hat(m) is a nonparametrically
#' calibrated version of the conditional mean. It uses the fitted values
#' from adaptive LASSO for the old data and new data, and uses a kernel
#' to smooth over them.
#' }
#'
#' \code{link_viol_sim} runs the simulation.
#'
#' @param nsims Number of simulations to run.
#' @param betas True effects, with first element equal to the intercept.
#'   Should be a length-(p+1) vector.
#' @param x_simulator Function whose first argument is n. Generates n
#'   replicates of X. The return value of this function should be an n x
#'   p matrix, or an n x 1 vector.
#' @param n Number of samples in each simulation.
#' @param true_link True link function. Should be an R function that
#'   takes one argument.
#' @param error_simulator Function whose first argument is n. Generates
#'   n replicates of epsilon. The return value of this function should
#'   be an n x 1 vector.
#' @param kern Kernel function to use for smoothing the predicted
#'   values.
#' @param testsize Sample size for the simulated validation dataset.
#'
#' @return \code{link_viol_sim} returns a list with two named elements:
#' \describe{
#' \item{betas}{nsims x p matrix of estimated coefficients for each
#' iteration of the simulation.}
#' \item{errors}{nsims x 3 matrix of three estimated errors.}
#' }
#'
#' \code{sim_data} returns a list with two named elements:
#' \describe{
#' \item{xs}{n x p matrix of predictors}
#' \item{ys}{length-n vector of outcomes}
#' }
link_viol_sim <- function(nsims, betas, x_simulator, n,
                          error_simulator = rnorm, kern = dnorm,
                          testsize = 5000) {

    stopifnot(length(betas) > 1)
    betas_est <- matrix(NA, nrow = nsims, ncol = length(betas))
    errors <- matrix(NA, nrow = nsims, ncol = 3)

    # only simulate the independent test dataset once
    test_dat <- sim_data(betas, x_simulator, error_simulator, testsize)
    for (i in 1:nsims) {
        # simulate data
        train_dat <- sim_data(betas, x_simulator, error_simulator, n)
        # fit adaptive lasso model
        alasso_cv <- adaptive_lasso(train_dat$xs, train_dat$ys)
        alasso_betas <- coef(alasso_cv, s = "lambda.min")

        # fit true model
        true_mod_betas <- fit_true_model(train_dat$xs, train_dat$ys)

        # get errors on independent test dataset
        true_preds <- true_mod_betas[1] +
            pnorm(test_dat$xs%*% true_mod_betas[-1])
        naive_preds <- predict(alasso_cv, newx = test_dat$xs,
                               s = "lambda.min")
        # calibration step
        calibrate_preds <- np_calibrate(test_dat$xs,
                                        train_dat$xs, train_dat$ys,
                                        alasso_cv, kern)
        # calculate errors
        errors[i, ] <- vapply(list(true_preds, naive_preds,
                                  calibrate_preds),
                             function(x) { mean((test_dat$ys - x)^2) },
                             FUN.VALUE = 2.1)
        betas_est[i, ] <- as.vector(alasso_betas)

    }
    return(list(betas = betas_est, errors = errors))
}

#' @describeIn link_viol_sim Simulate the data
sim_data <- function(betas, x_simulator, error_simulator, n) {
    xs <- scale(x_simulator(n), center = TRUE, scale = TRUE)
    stopifnot(dim(xs) == c(n, length(betas) - 1))
    errors <- error_simulator(n)
    ys <- betas[1] + pnorm(xs %*% betas[-1]) + errors
    return(list(xs = xs, ys = ys))
}

#' Nonparametric calibration step for predicting
#'
#' Nonparametric calibration step for predicting the mean
#'
#' The formula in the assignment does not include intercepts, but we do
#' include the intercepts here by using the \code{predict} function.
#' However, because the intercepts will be the same, after subtraction
#' they cancel out.
#'
#' @param new_xs Matrix of new predictors
#' @param xs Matrix of old predictors
#' @param ys Vector of old observed values
#' @param glmnet_obj A glmnet object for which there is a predict method
#' @param kern Kernel used for smoothing
#'
#' @return Nonparametrically calibrated predicted means.
np_calibrate <- function(new_xs, xs, ys, glmnet_obj, kern) {
    stopifnot(ncol(xs) == ncol(new_xs))
    yhat <- as.vector(predict(glmnet_obj, xs, s = "lambda.min"))
    yhat_new <- as.vector(predict(glmnet_obj, new_xs, s = "lambda.min"))
    kern_diffs <- kern(outer(yhat, yhat_new, `-`))
    mhat <- (t(kern_diffs) %*% ys) / colSums(kern_diffs)
    return(mhat)
}

#' Simulate correlated uniform random variates
#'
#' Simulate multiple iid correlated uniform random variates.
#'
#' First, with a given correlation matrix corr, we simulate
#'
#' (X_{i1}, ..., X_{ip}) ~ MVN(0, corr)
#'
#' for i = 1, 2, ... n. Then, we obtain
#'
#' (U_{i1}, ..., U_{ip}) = (F(X_{i1}), ..., F_(X_{ip})),
#'
#' where F is the standard normal PDF (note marginally, X_{i1} ~ N(0, 1)
#' since corr is assumed to have 1 on the diagonals). Then, the U will
#' marginally be uniformly distributed on \[0, 1\], with some
#' correlation.
#'
#' @param n Number of random variates to generate.
#' @param corr Correlation between the random variates. Must be a p x p
#' diagonal matrix, with 1 on the diagonal.
#'
#' @return n x p matrix of correlated uniform random variates, where
#' each row is iid.
runif_corr <- function(n, corr) {
    stopifnot(diag(corr) == rep(1, nrow(corr)))
    normals <- mvtnorm::rmvnorm(n, sigma = corr)
    unifs <- pnorm(normals)
    return(unifs)
}

#' Fit adaptive lasso model
#'
#' Fit adaptive lasso model with intercept. The penalty for each beta_j
#' coefficient in the model is weighted by 1 / |beta_ini_j|^gam, where
#' beta_ini_j is an initial estimate of beta, here obtained through OLS.
#' Then, a LASSO is run. The optimal lambda is selected through
#' cross-validation.
#'
#' @param xs Matrix of predictors. Should be standardized and centered
#'   to have mean 0, variance 1
#' @param ys Matrix or vector of observations. Does not need to be
#'   centered.
#' @param gam Power for penalty weights.
#' @param nfolds Number of cross-validation folds.
#'
#' @return A cv.glmnet object
#'
adaptive_lasso <- function(xs, ys, gam = 1, nfolds = 10) {
    ols_reg <- lm(ys ~ xs)
    beta_ini <- coef(ols_reg)[-1]
    penalty_weights = 1 / (abs(beta_ini)^gam)
    cv <- glmnet::cv.glmnet(xs, ys, standardize = FALSE,
                            intercept = TRUE, nfolds = nfolds,
                            penalty.factor = penalty_weights)
    return(cv)
}
