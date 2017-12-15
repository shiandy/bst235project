library(bst235Project)
library(purrr)
library(mvtnorm)
library(tools)
library(foreach)
library(parallel)
library(doParallel)

# HELPER FUNCTIONS -----------------------------------------------------

# Generate an exchangeable correlation matrix. Diagonal elements are 1,
# off-diagonal elements are rho
exch_corr <- function(n, rho) {
    stopifnot(rho < 1)
    ret <- matrix(rho, nrow = n, ncol = n)
    diag(ret) <- 1
    return(ret)
}

# Generate block exchangeable correlation matrix. n1 and n2 are the
# sizes of the two blocks, so the resultant correlation matrix is (n1 +
# n2) x (n1 + n2). rho1 is the pairwise correlation in block 1, rho2 is
# the pairwise correlation in block 2, and rho12 is the correlation
# between elements in block 1 and block 2.
block_corr <- function(n1, n2, rho1, rho12, rho2) {
    stopifnot(all(c(rho1, rho12, rho2) < 1))
    block1 <- exch_corr(n1, rho1)
    block2 <- exch_corr(n2, rho2)
    block12 <- matrix(rho12, nrow = n1, ncol = n2)
    ret <- rbind(cbind(block1, block12), cbind(t(block12), block2))
    return(ret)
}

rmvnorm_distorted <- function(n, corr_mat) {
    xs <- mvtnorm::rmvnorm(n, sigma = corr_mat)
    xs[,8] <- xs[,8] + 0.5*(xs[,7]^2-1)
    return(xs)
}

# SIMULATION SETUP -----------------------------------------------------

# number of simulations to run
nsims <- 1000
ns <- c(500, 1000, 5000)
rho1s <- c(0, 0.3)
rho2s <- c(0, 0.3)
x_dist <- c("normal", "uniform", "distorted_normal")

n_zeros <- 5
n_nonzeros <- 4
betas <- c(3, rep(0, n_zeros), -0.3, -0.1, 0.1, 0.3)
sds <- c(0.01, 0.3)

sim_params <- expand.grid(rho1 = rho1s, rho2 = rho2s, x_dist = x_dist,
                          n = ns, sd = sds)
sim_params$rho12 <-
    ifelse(sim_params$rho1 == 0.3 & sim_params$rho2 == 0.3, 0.3, 0)

# RUN THE SIMULATION ---------------------------------------------------

# file to save the results
data_csv <- "../data/simdata.csv"

# parallelize
registerDoParallel(cores = 4)
#res_df_lst <- vector("list", nrow(sim_params))
system.time({
    res_df <- foreach(i = 1:nrow(sim_params),
                      .packages = c("bst235Project", "tools", "purrr"),
                      .combine = rbind) %dopar% {
        set.seed(1263)
        # extract parameters
        n <- sim_params$n[i]
        rho1 <- sim_params$rho1[i]
        rho12 <- sim_params$rho12[i]
        rho2 <- sim_params$rho2[i]
        x_dist <- sim_params$x_dist[i]
        sd <- sim_params$sd[i]

        error_simulator <- function(n) {
            return(rnorm(n, sd = sd))
        }

        print(sprintf("%s... running iteration %d of %d; n = %d",
              Sys.time(), i, nrow(sim_params), n))
        corr_mat <- block_corr(n_zeros, n_nonzeros, rho1, rho12, rho2)
        if (x_dist == "normal") {
            x_simulator <- partial(rmvnorm, sigma = corr_mat)
        } else if (x_dist == "uniform") {
            x_simulator <- partial(runif_corr, corr = corr_mat)
        }
        else if (x_dist == "distorted_normal") {
            x_simulator <- partial(rmvnorm_distorted,
                                   corr_mat = corr_mat)
        }
        sim_ret <- link_viol_sim(nsims, betas, x_simulator, n,
                                 error_simulator, testsize = 10000,
                                 cv = FALSE)

        # get the result, add data about simulation parameters
        df_cur <- as.data.frame(cbind(sim_ret$betas, sim_ret$errors))
        df_cur$n <- n
        df_cur$rho1 <- rho1
        df_cur$rho12 <- rho12
        df_cur$rho2 <- rho2
        df_cur$x_dist <- x_dist
        df_cur$sd <- sd
        #res_df_lst[[i]] <- df_cur
        df_cur
    }
})
#res_df <- do.call(rbind, res_df_lst)
write.csv(res_df, file = data_csv, row.names = FALSE)
print(sprintf("Simulation data written to %s, with md5sum of: %s",
              data_csv, md5sum(data_csv)))
