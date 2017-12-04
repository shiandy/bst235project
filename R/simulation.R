#' Simulate violations of link function
#'
#' Extra description.
#'
#' @param nsims Number of simulations to run.
#' @param beta True effects. Should be a length-p vector.
#' @param x_simulator Function whose first argument is n. Generates n
#'   replicates of X. The return value of this function should be an n x
#'   p matrix, or an n x 1 vector.
#' @param true_link True link function. Should be an R function that
#'   takes one argument.
#' @param wrong_link Misspecified link function. Should be an R function
#'   that takes one argument.
#' @param error_simulator Function whose first argument is n. Generates
#'   n replicates of epsilon. The return value of this function should
#'   be an n x 1 vector.
#' @param ns Vector of n to use.
#' @param lambdas Vector of lambdas to use. Should be the same length as
#'   ns.
#' @param kern Kernel function to use for smoothing the predicted
#'   values.
#' @param bandwidths Vector of bandwidths to use for smoothing the
#'   predicted values.
#'
#' @return A list with two named elements
#' \describe{
#' \item{betas}{nsims x p DataFrame of estimated coefficients for each
#' iteration of the simulation.}
#' \item{errors}{nsims x 3 DataFrame of three estimated errors.}
#' }
#'
#' @examples
#' 1 + 1
link_viol_sim <- function(nsims, betas, x_simulator, true_link,
                          error_simulator = rnorm, n, kern = dnorm) {
    # simulate data
    train_dat <- sim_data(betas, x_simulator, true_link,
                          error_simulator, n) 
    # fit adaptive lasso model
    alasso_cv <- adaptive_lasso(train_dat$xs, train_dat$ys)
    alasso_betas <- coef(alasso_cv, s = "lambda.min")

    # calibration step
    test_dat <- sim_data(betas, x_simulator, true_link,
                          error_simulator, n)
    # fit true model
    glm(train_dat$ys ~ train_dat$xs,
        family = gaussian(link = make.link("probit")))
    glm(train_dat$ys ~ train_dat$xs,
        family = gaussian())

    # calculate errors
}

sim_data <- function(betas, x_simulator, true_link,
                     error_simulator, n) {
    xs <- scale(x_simulator(n), center = TRUE, scale = TRUE)
    stopifnot(dim(xs) == c(n, length(betas)))
    errors <- error_simulator(n)
    ys <- scale(true_link(xs %*% betas) + errors,
                center = TRUE, scale = FALSE)
    return(list(xs = xs, ys = ys))
}


x_mvnorm <- function(n) {
    return(rmvnorm(n, sigma = diag(2)))
}

np_calibrate <- function(new_xs, new_ys, xs, ys, glmnet_obj, kern) {
    yhat <- as.vector(predict(glmnet_obj, xs, s = "lambda.min"))
    yhat_new <- as.vector(predict(glmnet_obj, new_xs, s = "lambda.min"))
    kern_diffs <- kern(outer(yhat, yhat_new, `-`))
    mhat <- (t(kern_diffs) %*% ys) / colSums(kern_diffs)

}

#' Simulate correlated uniform random variates
#'
#' Simulate multiple iid correlated uniform random variates.
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
#' Fit adaptive lasso model
#'
#' @param xs Matrix of predictors. Should be standardized and centered
#' to have
#' @param ys Matrix or vector of observations.
adaptive_lasso <- function(xs, ys, gam = 1, nfolds = 10) {
    ols_reg <- lm(ys ~ xs - 1)
    beta_ini <- coef(ols_reg)
    penalty_weights = 1 / (abs(beta_ini)^gam)
    cv <- glmnet::cv.glmnet(xs, ys, standardize = FALSE,
                            intercept = FALSE, nfolds = nfolds,
                            penalty.factor = penalty_weights)
    return(cv)
}
