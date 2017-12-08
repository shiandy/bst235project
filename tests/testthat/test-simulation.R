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
