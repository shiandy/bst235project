context("test-simulation.R")

# setup code
betas <- c(1, 0, 5)
x_simulator <- function(n) {
    mvtnorm::rmvnorm(n, sigma = diag(2))
}
error_simulator <- function(n) { rnorm(n, sd = 0.5) }
nsims <- 10
n <- 100
set.seed(1262)
system.time({
    sim_res <- link_viol_sim(nsims, betas, x_simulator, n,
                             error_simulator = error_simulator)
})

test_that("link_viol_sim returns objects with correct dimensions", {
    expect_equal(dim(sim_res$betas), c(nsims, length(betas)))
    expect_equal(dim(sim_res$errors), c(nsims, 3))
})
