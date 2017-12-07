context("test-simulation.R")

# setup code
n_zeros <- 5
n_nonzeros <- 3
#betas <- c(3, rep(0, n_zeros), 0.5, 1, 1.5)
betas <- c(3, rep(0, n_zeros), -0.3, 0.1, 0.3)
corr_mat <- matrix(0.3, nrow = n_zeros + n_nonzeros,
                   ncol = n_zeros + n_nonzeros)
diag(corr_mat) <- 1
x_simulator <- function(n) {
    return(runif_corr(n, corr_mat))
}
error_simulator <- function(n) { rnorm(n, sd = 1) }
nsims <- 5
n <- 1000
set.seed(1262)
system.time({
    sim_res <- link_viol_sim(nsims, betas, x_simulator, n,
                             error_simulator = error_simulator,
                             testsize = 10000)
})

test_that("link_viol_sim returns objects with correct dimensions", {
    expect_equal(dim(sim_res$betas), c(nsims, length(betas)))
    expect_equal(dim(sim_res$errors), c(nsims, 3))
})

test_that("runif_corr has the correct bound", {
    unifs <- runif_corr(1000, corr_mat, min = -2, max = 2)
    expect_gte(min(unifs), -2)
    expect_lte(max(unifs), 2)
})


# debugging simulation
set.seed(1262)
test_dat <- sim_data(betas, x_simulator, error_simulator,
                     n = 10000)
train_dat <- sim_data(betas, x_simulator, error_simulator, n)
alasso_cv <- adaptive_lasso(train_dat$xs, train_dat$ys)
alasso_betas <- coef(alasso_cv, s = "lambda.min")

# fit true model
true_mod_betas <- fit_true_model(train_dat$xs, train_dat$ys)

true_preds <- true_mod_betas[1] +
    pnorm(test_dat$xs%*% true_mod_betas[-1])
naive_preds <- predict(alasso_cv, newx = test_dat$xs, s = "lambda.min")
calibrate_preds <- np_calibrate(test_dat$xs, train_dat$xs, train_dat$ys,
                                alasso_cv)
par(mfrow = c(2, 2))
hist(true_preds)
plot(true_preds, naive_preds)
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(true_preds, calibrate_preds)
abline(a = 0, b = 1, col = "red", lwd = 2)
hist(calibrate_preds)

