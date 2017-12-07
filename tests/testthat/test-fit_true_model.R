context("test-fit_true_model.R")

# setup code
betas <- c(1, 0, 5)
x_simulator <- function(n) {
    mvtnorm::rmvnorm(n, sigma = diag(1, nrow = 2, ncol = 2))
}
n <- 100

# simulate data with no errors
error_simulator0 <- function(n) { 0 }
dat0 <- sim_data(betas, x_simulator, error_simulator0, n)

test_that("Objective function is min at correct parameters", {
    expect_lte(objective(betas, dat0$xs, dat0$ys),
               objective(betas * 2, dat0$xs, dat0$ys))
    expect_lte(objective(betas, dat0$xs, dat0$ys),
               objective(betas - 2, dat0$xs, dat0$ys))
})

test_that("Gradient is 0 at correct parameters for no error data", {
    expect_equal(gradient(betas, dat0$xs, dat0$ys),
                 rep(0, length(betas)))
    grad <- gradient(betas * 2, dat0$xs, dat0$ys)
    expect_false(isTRUE(all.equal(grad,  rep(0, length(betas)))))
})

test_that("We recover the true parameters when there is no error", {
    fit_betas <- fit_true_model(dat0$xs, dat0$ys)
    expect_equal(fit_betas, betas, check.attributes = FALSE)
})
