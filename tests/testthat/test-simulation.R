context("test-simulation.R")

# setup code
n_zeros <- 5
n_nonzeros <- 3
#betas <- c(3, rep(0, n_zeros), 0.5, 1, 1.5)
betas <- c(3, rep(0, n_zeros), -0.3, 0.1, 0.5)
corr_mat <- matrix(0.3, nrow = n_zeros + n_nonzeros,
                   ncol = n_zeros + n_nonzeros)
diag(corr_mat) <- 1
x_simulator <- function(n) {
    return(runif_corr(n, corr_mat))
}
error_simulator <- function(n) { rnorm(n, sd = 0.3) }
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
x_simulator2 <- function(n) {
    return(mvtnorm::rmvnorm(n, sigma = corr_mat))
}
set.seed(1262)
test_dat <- sim_data(betas, x_simulator2, error_simulator,
                     n = 10000)
train_dat <- sim_data(betas, x_simulator2, error_simulator, n)
alasso_cv <- adaptive_lasso(train_dat$xs, train_dat$ys)
alasso_betas <- coef(alasso_cv, s = "lambda.min")

# fit true model
true_mod_betas <- fit_true_model(train_dat$xs, train_dat$ys)

# get predictions
true_preds <- true_mod_betas[1] +
    pnorm(test_dat$xs%*% true_mod_betas[-1])
naive_preds <- predict(alasso_cv, newx = test_dat$xs, s = "lambda.min")

# do calibration predictions
yhat <- predict(alasso_cv, newx = train_dat$xs, s = "lambda.min")
best_h <- select_bandwidth(yhat, train_dat$ys, nfolds = 10)
calibrate_preds <- np_calibrate(test_dat$xs, train_dat$xs, train_dat$ys,
                                alasso_cv, bandwidth = best_h)
errors_cur <- vapply(list(true_preds, naive_preds, calibrate_preds),
                     function(x) { mean((test_dat$ys - x)^2,
                                        na.rm = TRUE) }, FUN.VALUE = 2.1)
print(errors_cur)

# histogram of predicted and true values
par(mfrow = c(2, 2))
hist(test_dat$ys)
hist(true_preds)
hist(naive_preds)
hist(calibrate_preds)

# scatterplots
plt_lims <- c(min(test_dat$ys), max(test_dat$ys))
par(mfrow = c(2, 2))
plot(test_dat$ys, true_preds, ylim = plt_lims)
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(test_dat$ys, naive_preds, ylim = plt_lims)
abline(a = 0, b = 1, col = "red", lwd = 2)
plot(test_dat$ys, calibrate_preds, ylim = plt_lims)
abline(a = 0, b = 1, col = "red", lwd = 2)
hist(test_dat$xs %*% true_mod_betas[-1])

test_that("Cross validation", {
    skip("Run this yourself... takes too long")
    system.time({
        set.seed(1262)
        sim_res_cv <- link_viol_sim(1000, betas, x_simulator, n,
                                     error_simulator = error_simulator,
                                     testsize = 10000, cv = TRUE)
    })
    system.time({
        set.seed(1262)
        sim_res_nocv <- link_viol_sim(1000, betas, x_simulator, n,
                                     error_simulator = error_simulator,
                                     testsize = 10000, cv = FALSE)
    })
    boxplot(sim_res_cv$errors)
    boxplot(sim_res_nocv$errors)
    apply(sim_res_cv$errors, 2, median)
    apply(sim_res_nocv$errors, 2, median)
})
