library(bst235Project)
library(tools)
library(foreach)
library(parallel)
library(doParallel)

########################################################################
# HELPER FUNCTIONS
########################################################################

# Generate an exchangeable correlation matrix.
exch_corr <- function(n, rho) {
    stopifnot(rho < 1)
    ret <- matrix(rho, nrow = n, ncol = n)
    diag(ret) <- 1
    return(ret)
}

# generate block correlation matrix
block_corr <- function(n1, n2, rho1, rho12, rho2) {
    stopifnot(all(c(rho1, rho12, rho2) < 1))
    block1 <- exch_corr(n1, rho1)
    block2 <- exch_corr(n2, rho2)
    block12 <- matrix(rho12, nrow = n1, ncol = n2)
    ret <- rbind(cbind(block1, block12), cbind(t(block12), block2))
    return(ret)
}

# generate x_simulator for X ~ multivariate normal
rmvnorm_generator <- function(corr_mat) {
    ret <- function(n) {
        return(mvtnorm::rmvnorm(n, sigma = corr_mat))
    }
    return(ret)
}

# generate x_simulator for X ~ correlated uniform
runif_corr_generator <- function(corr_mat) {
    ret <- function(n) {
        return(runif_corr(n, corr_mat))
    }
    return(ret)
}

########################################################################
# SIMULATION SETUP
########################################################################

# number of simulations to run
nsims <- 2
ns <- c(1000, 10000)
rho1s <- c(0, 0.3)
rho2s <- c(0, 0.3)
x_dist <- c("normal", "uniform")

n_zeros <- 5
n_nonzeros <- 3
betas <- c(3, rep(0, n_zeros), 1, 3, 2)

error_simulator <- function(n) {
    return(rnorm(n, sd = 1))
}

sim_params <- expand.grid(rho1 = rho1s, rho2 = rho2s, x_dist = x_dist,
                          n = ns)
sim_params$rho12 <-
    ifelse(sim_params$rho1 == 0.3 & sim_params$rho2 == 0.3, 0.3, 0)

########################################################################
# RUN THE SIMULATION
########################################################################

# file to save the results
data_csv <- "../data/simdata.csv"

registerDoParallel(cores = 4)
#res_df_lst <- vector("list", nrow(sim_params))
system.time({
    res_df <- foreach(i = 1:nrow(sim_params), .combine = rbind,
                      .packages = c("bst235Project", "tools")) %dopar% {
        set.seed(1262)
        # extract parameters
        n <- sim_params$n[i]
        rho1 <- sim_params$rho1[i]
        rho12 <- sim_params$rho12[i]
        rho2 <- sim_params$rho2[i]
        x_dist <- sim_params$x_dist[i]

        print(sprintf("%s... running iteration %d of %d; n = %d",
              Sys.time(), i, nrow(sim_params), n))
        corr_mat <- block_corr(n_zeros, n_nonzeros, rho1, rho12,
                               rho2)
        if (x_dist == "normal") {
            x_simulator <- rmvnorm_generator(corr_mat)
        } else if (x_dist == "uniform") {
            x_simulator <- runif_corr_generator(corr_mat)
        }
        sim_ret <- link_viol_sim(nsims, betas, x_simulator, n,
                                 error_simulator, testsize = 10000)

        df_cur <- as.data.frame(cbind(sim_ret$betas,
                                      sim_ret$errors))
        df_cur$n <- n
        df_cur$rho1 <- rho1
        df_cur$rho12 <- rho12
        df_cur$rho2 <- rho2
        df_cur$x_dist <- x_dist
        #res_df_lst[[i]] <- df_cur
        df_cur
    }
})
#res_df <- do.call(rbind, res_df_lst)
write.csv(res_df, file = data_csv, row.names = FALSE)
print(sprintf("Simulation data written to %s, with md5sum of: %s",
              data_csv, md5sum(data_csv)))
