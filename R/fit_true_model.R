#' Fit true regression model
#'
#' Fit true regression model E(Y_i | X_i) = alpha_0 + Phi(beta^T X_i),
#' where Phi is the standard Normal CDF.
#'
#' Under the hood, the fitting is done using BFGS, as implemented in the
#' \code{optim} function in R.
#'
#' @param xs n x p matrix of predictors
#' @param ys Length-n vector of outcomes
#' @param init Initial guess for parameters
#'
#' @return \code{fit_true_model} returns a length-(p + 1) vector of
#'   coefficients, where the first is the intercept and the rest are the
#'   beta values. If the fitting algorithm did not converge,
#'   \code{fit_true_model} will fail with an error.
#'
#'   \code{objective} and \code{gradient} return the value of the
#'   objective function and gradient evaluated for some coefficients
#'   alpha_beta, respectively.
fit_true_model <- function(xs, ys, init = NULL) {
    stopifnot(nrow(xs) == nrow(ys))
    #nparams <- ncol(xs) + 1
    if (is.null(init)) {
        init <- coef(lm(ys ~ xs))
    }
    solution <- optim(init, fn = objective, gr = gradient,
                      method = "BFGS", xs = xs, ys = ys,
                      control = list(maxit = 10000))
    if (solution$convergence) {
        stop("Fitting algorithm for true model failed to converge.")
    }
    return(solution$par)
}

#' @param alpha_betas Value of coefficients to evaluate the function.
#'   First entry must be the intercept, and the rest are the other
#'   coefficients.
#' @describeIn fit_true_model Evaluate objective function
objective <- function(alpha_betas, xs, ys) {
    stopifnot(nrow(xs) == nrow(ys))
    stopifnot(ncol(xs) == length(alpha_betas) - 1)
    alpha <- alpha_betas[1]
    betas <- alpha_betas[-1]
    ret <- 0.5 * sum((ys - alpha - pnorm(xs %*% betas))^2)
    return(ret)
}

#' @describeIn fit_true_model Evaluate objective function's gradient
gradient <- function(alpha_betas, xs, ys) {
    stopifnot(nrow(xs) == nrow(ys))
    stopifnot(ncol(xs) == length(alpha_betas) - 1)
    alpha <- alpha_betas[1]
    betas <- alpha_betas[-1]
    grad <- rep(NA, length(alpha_betas))
    linear_term <- xs %*% betas
    sum_term <- ys - alpha - pnorm(linear_term)
    grad[1] <- -sum(sum_term)
    grad[2:length(alpha_betas)] <-
        -t(xs) %*% (sum_term * dnorm(linear_term))
    return(grad)
}
