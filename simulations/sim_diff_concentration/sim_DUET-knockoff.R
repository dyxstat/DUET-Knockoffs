


# # install.packages("devtools")
# devtools::install_github("dyxstat/DUET-Knockoffs")

library(parallel)
library(MASS)
library(gtools)
library(dirmult)
library(writexl)
library(readxl)
library(HMP)

Sys.setenv(
  MC_CORES = "1",           
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

RNGkind("Mersenne-Twister")
set.seed(42)

B <- 100     # 100 times monte carlo

n_confounders <- 3

filter.thresh <- 0.2
n_rej_stop <- 100

have_bias <- 1; have_diff_bias <- 1; 
bias_sd1 <- 0.1; bias_sd2 <- 0.1

depth_fold <- 10
depth.mu1 <- 100000; depth.mu2 <- 100000 * depth_fold

depth.sd1 <- depth.mu1/5; depth.sd2 <- depth.mu2/5 
depth.lower <- 10000

df <- readRDS("ZhuF_count.rds")  # YeZ_count, Dha_count
df <- df[, colSums(df) > 0]

pi_est <- DM.MoM(df)$pi
n.otus <- length(pi_est)

## -----------------------------
## 100 times Monte Carlo over grid (DUET-knockoff)
## -----------------------------

c_grid <- c(1e4, 1e5, 1e6)
beta_grid <- c(4, 6, 8, 10, 12, 14)
prop.diff_grid <- c(0.15)
nsam_grid <- list(c(50, 50)) 
fdr_grid <- c(0.2)
## -----------------------------
## -----------------------------
## -----------------------------


############# estimate phi #######################

gmean <- function(x) exp(mean(log(x[x > 0])))  ## GM(d)

size_factors_from_depth <- function(depth) depth / gmean(depth)  ## si

phi_from_vector <- function(y, depth, trim = 0.1, eps = 1e-8) {
  s  <- size_factors_from_depth(depth)
  m  <- sum(y) / sum(s)
  mu <- pmax(s * m, eps)
  w  <- ((y - mu)^2 - mu) / (mu^2)
  w  <- w[is.finite(w)]
  if (!length(w)) return(list(phi=0, mu=mu, s=s, m_hat=m))
  w <- sort(w)
  k <- floor(length(w) * trim)
  if (k > 0 && length(w) > 2*k) w <- w[(k+1):(length(w)-k)]
  list(phi = max(mean(w), 0), mu = mu, s = s, m_hat = m)
}

# nb_zero_prob <- function(mu, phi) (1 + phi * mu)^(-1/phi)  ## nb zero prob

depth <- rowSums(df, na.rm = TRUE)
s_vec <- size_factors_from_depth(depth)
sum_s <- sum(s_vec)
trim_phi <- 0.1; eps_phi <- 1e-8
phi_cores <- max(1L, parallel::detectCores() - 1L)

phi <- unlist(parallel::mclapply(
  X = seq_len(n.otus),
  FUN = function(i) {
    y  <- df[, i]
    m  <- sum(y) / sum_s
    mu <- pmax(s_vec * m, eps_phi)
    w  <- ((y - mu)^2 - mu) / (mu^2)
    w  <- w[is.finite(w)]
    if (!length(w)) return(0)
    w <- sort(w)
    k <- floor(length(w) * trim_phi)
    if (k > 0 && length(w) > 2*k) w <- w[(k+1):(length(w)-k)]
    max(mean(w), 0)
  },
  mc.cores = phi_cores,
  mc.preschedule = TRUE
))

##########################################################


## -----------------------------
## single simulation (c, beta, prop.diff, nsam, fdr_target, seed)
## return FDR and Power per method
## -----------------------------
simulate_one_param <- function(c_val, beta_val, prop.diff_val, nsam_val, fdr_target, seed) {
  
  n_cores <- 1
  
  # ==== Begin of Data Preparation ==== #
  
  n_DA <- ceiling(n.otus * prop.diff_val)
  
  n.sam.grp1 <- nsam_val[1]
  n.sam.grp2 <- nsam_val[2]
  
  n_sam <- n.sam.grp1 + n.sam.grp2
  
  ## 1) random choose DA otus
  causal.otus.idx <- sample(which(pi_est > 0), n_DA)
  noncausal.otus.idx <- setdiff(1:n.otus, causal.otus.idx)
  
  beta.otu <- rep(beta_val, length(causal.otus.idx))
  beta.otu.log <- log(beta.otu)
  
  # 2)
  if (n_confounders > 0) {
    
    betaC <- log(1.1)
    
    confounding.otus <- matrix(NA, nrow=n_confounders, ncol=5)
    w <- which(pi_est >= 0.005)
    for (r in 1:n_confounders) confounding.otus[r,] <- sort(sample(w, 5))
  }
  
  # 3) define bias factor
  if (have_bias) {
    
    bias.factor1.log <- rnorm(n.otus, 0, bias_sd1)
    bias.factor2.log <- rnorm(n.otus, 0, bias_sd2)
    
    if (have_diff_bias) {
      sub1 <- 1:5
      sub2 <- 11:15
      
      top5 <- order(pi_est, decreasing = TRUE)[1:5]
      sub1 <- c(causal.otus.idx[sub1], 
                setdiff(sample(noncausal.otus.idx, length(noncausal.otus.idx)/5), top5))
      sub2 <- c(causal.otus.idx[sub2], 
                setdiff(sample(noncausal.otus.idx, length(noncausal.otus.idx)/5), top5))
      
      bias.factor1.log[sub1] <- -5
      bias.factor2.log[sub2] <- -5
    }
    bias.factor1 <- exp(bias.factor1.log)
    bias.factor2 <- exp(bias.factor2.log)
    top_taxa <- order(pi_est, decreasing = TRUE)[1]
    bias.factor1[top_taxa] <- 1
    bias.factor2[top_taxa] <- 1
  }
  
  ## 4) generate Covariate C, simulated library size
  if (n_confounders > 0) {
    C <- matrix(NA, nrow = n_sam, ncol = n_confounders)
    for (r in 1:n_confounders) {
      C[, r] <- c(runif(n.sam.grp1, -1, 1), runif(n.sam.grp2, 0, 2))
    }
  } else C <- NULL
  
  # transform moments of normal distirbution to lognormal
  cv1 <- depth.sd1 / depth.mu1
  cv2 <- depth.sd2 / depth.mu2
  
  meanlog1 <- log(depth.mu1) - 0.5 * log(1 + cv1^2)
  sdlog1 <- sqrt(log(1 + cv1^2))
  
  meanlog2 <- log(depth.mu2) - 0.5 * log(1 + cv2^2)
  sdlog2 <- sqrt(log(1 + cv2^2))
  
  depth1.sim <- rlnorm(n_sam, meanlog1, sdlog1)
  depth1.sim[depth1.sim < depth.lower] <- depth.lower  # trunct
  depth1.sim <- round(depth1.sim)  # rounding
  
  depth2.sim <- rlnorm(n_sam, meanlog2, sdlog2)
  depth2.sim[depth2.sim < depth.lower] <- depth.lower  # trunct
  depth2.sim <- round(depth2.sim)  # rounding
  
  ## 5) Define binary trait Y, baseline relative abundance matrix pi.table.sim
  Y <- c(rep(0, n.sam.grp1), rep(1, n.sam.grp2)) 
  
  # Dirichlet sampling
  
  c <- c_val          # 1e4, 1e5, 1e6
  alpha <- c * pi_est
  
  pi.table.sim <- rdirichlet(n = n_sam, alpha)
  
  rownames(pi.table.sim) <- paste0("sub", 1:n_sam)
  colnames(pi.table.sim) <- paste0("taxon", 1:n.otus)
  
  pi.table1.sim <- pi.table.sim  # 16s
  pi.table2.sim <- pi.table.sim  # shotgun
  
  causal.otus <- colnames(pi.table.sim)[causal.otus.idx]          ## causal names
  noncausal.otus <- setdiff(colnames(pi.table.sim), causal.otus)  ## noncausal names
  
  ################################################################
  # 6) Introduce causal effect
  pi.table1.sim[, causal.otus.idx] <- pi.table1.sim[, causal.otus.idx] * exp(Y %*% t(beta.otu.log))
  pi.table2.sim[, causal.otus.idx] <- pi.table2.sim[, causal.otus.idx] * exp(Y %*% t(beta.otu.log))
  
  # Introduce confounding effect
  if (n_confounders > 0) {
    for (r in 1:n_confounders) {
      # betaC^C[, r] = exp(C * log(betaC))
      conf_factor <- exp(betaC * C[, r])    # or  (1.1)^C[, r]
      pi.table1.sim[, confounding.otus[r,]] <- sweep(
        pi.table1.sim[, confounding.otus[r,], drop = FALSE],
        1, conf_factor, "*"
      )
      pi.table2.sim[, confounding.otus[r,]] <- sweep(
        pi.table2.sim[, confounding.otus[r,], drop = FALSE],
        1, conf_factor, "*"
      )
    }
  }
  
  # Introduce bias effect
  if (have_bias == 1) {
    pi.table1.sim <- sweep(pi.table1.sim, 2, bias.factor1, "*")
    pi.table2.sim <- sweep(pi.table2.sim, 2, bias.factor2, "*")
  }
  
  # # Normalization
  # pi.table1.sim <- pi.table1.sim / rowSums(pi.table1.sim)
  # pi.table2.sim <- pi.table2.sim / rowSums(pi.table2.sim)
  
  ################################################################
  
  # 7) Poisson-Gamma sampling
  mu1 <- pi.table1.sim * depth1.sim    # mean absolute abundance matrix for 16s
  mu2 <- pi.table2.sim * depth2.sim
  
  
  phi1_safe <- pmax(phi, 1e-8)
  phi2_safe <- pmax(phi, 1e-8)
  
  shape1 <- 1 / phi1_safe
  shape2 <- 1 / phi2_safe
  
  otu.table1.sim <- matrix(0, n_sam, n.otus)
  otu.table2.sim <- matrix(0, n_sam, n.otus)
  
  idx <- (phi >= 1e-8)     # Poisson-gamma for phi != 0
  idx_p  <- !idx           # Poisson sampling for phi = 0
  
  # 16S count matrix
  if (any(idx)) {
    phi_pg1 <- phi1_safe[idx]
    mu1_pg <- mu1[, idx, drop = FALSE]
    A1 <- rgamma(n = n_sam * sum(idx),
                 shape = rep(1 / phi_pg1, each = n_sam),
                 rate = rep(1 / phi_pg1, each = n_sam) / c(mu1_pg))
    X1 <- rpois(n_sam * sum(idx), lambda = A1)
    otu.table1.sim[, idx] <- matrix(X1, nrow = n_sam)
  }
  if (any(idx_p)) {
    mu1_p <- mu1[, idx_p, drop = FALSE]
    x_p_1 <- rpois(n = n_sam * sum(idx_p), lambda = c(mu1_p))
    otu.table1.sim[, idx_p] <- matrix(x_p_1, nrow = n_sam)
  }
  
  # SMS count matrix
  if (any(idx)) {
    phi_pg2 <- phi2_safe[idx]
    mu2_pg <- mu2[, idx, drop = FALSE]
    A2 <- rgamma(n = n_sam * sum(idx),
                 shape = rep(1 / phi_pg2, each = n_sam),
                 rate = rep(1 / phi_pg2, each = n_sam) / c(mu2_pg))
    X2 <- rpois(n_sam * sum(idx), lambda = A2)
    otu.table2.sim[, idx] <- matrix(X2, nrow = n_sam)
  }
  if (any(idx_p)) {
    mu2_p <- mu2[, idx_p, drop = FALSE]
    x_p_2 <- rpois(n = n_sam * sum(idx_p), lambda = c(mu2_p))
    otu.table2.sim[, idx_p] <- matrix(x_p_2, nrow = n_sam)
  }
  
  colnames(otu.table1.sim) <- colnames(otu.table2.sim) <- paste0("taxon", 1:n.otus)
  rn <- sprintf("id_%03d", seq_len(nrow(otu.table1.sim)))
  rownames(otu.table1.sim) <- rn
  rownames(otu.table2.sim) <- rn
  
  ## 8) Filter & pool
  prop.presence1 <- colMeans(otu.table1.sim > 0)
  otus.keep1 <- which(prop.presence1 >= filter.thresh)
  otu.table1.sim.filter <- otu.table1.sim[, otus.keep1, drop=FALSE]
  
  prop.presence2 <- colMeans(otu.table2.sim > 0)
  otus.keep2 <- which(prop.presence2 >= filter.thresh)
  otu.table2.sim.filter <- otu.table2.sim[, otus.keep2, drop=FALSE]
  
  otu.table.sim.pool <- otu.table1.sim + otu.table2.sim
  
  prop.presence.pool <- colMeans(otu.table.sim.pool > 0)
  otus.keep.pool <- which(prop.presence.pool >= filter.thresh)
  otu.table.sim.pool.filter <- otu.table.sim.pool[, otus.keep.pool, drop=FALSE]
  
  # ==== End of Data Preparation ==== #
  
  # ==== ZINBSK ==== #
  W1 <- otu.table1.sim[, otus.keep.pool, drop = FALSE]
  W2 <- otu.table2.sim[, otus.keep.pool, drop = FALSE]
  
  M1 <- rowSums(W1)
  M2 <- rowSums(W2)
  
  
  # data_x
  
  data_x_1 <- data_x_2 <- cbind(C[,1], C[,2], C[,3])
  data_x_1 <- as.data.frame(data_x_1)             # 16s data_x
  data_x_2 <- as.data.frame(data_x_2)             # SMS data_x
  colnames(data_x_1) <- c("X_1", "X_2", "X_3")
  colnames(data_x_2) <- c("X_1", "X_2", "X_3")
  
  ##
  
  count_16s <- cbind(data_x_1, M = M1, W1)
  count_SMS <- cbind(data_x_2, M = M2, W2)
  count_SK <- rbind(count_16s, count_SMS)
  
  # define parameter
  n_group <- 2
  n_data <- 2
  class_K <- factor(rep(1:n_data, each = n_sam))  # label for data source
  
  y <- rep(c(rep(1, n.sam.grp1), rep(2, n.sam.grp2)), times = n_data) # y
  
  # select column
  n_x1 <- 3
  data_x <- as.data.frame(count_SK[, c(1:n_x1)])
  W <- as.data.frame(count_SK[, -c(1:(n_x1 + 1))])
  M <- count_SK[, n_x1 + 1]
  
  T_var <- which(colnames(W) %in% causal.otus)
  
  res.DUETknockoff <- DUET_knockoff(
    W = W, class_K = class_K, data_x = data_x, M = M, y = y, T_var = T_var,
    fdr = fdr_target, test_statistic = "DE", filter_statistics = 3,
    offset=1
  )
  
  duetknockoff.fdr <- res.DUETknockoff$FDRPower[2]
  duetknockoff.power <- res.DUETknockoff$FDRPower[3]
  
  fdr_vec <- DUETknockoff=duetknockoff.fdr
  
  power_vec <- DUETknockoff=duetknockoff.power
  
  list(fdr = fdr_vec, power = power_vec)
}


## store combination K
Kc <- length(c_grid); Kb <- length(beta_grid); Kf <- length(fdr_grid); 
Kp <- length(prop.diff_grid); Kn <- length(nsam_grid)

K <- Kc * Kb * Kf * Kp * Kn

## combo indenx and labels
idx5 <- function(ic, ib, ip, i_nsam, ifdr) {
  ((((ic - 1) * Kb + (ib - 1)) * Kp + (ip - 1)) * Kn + (i_nsam - 1)) * Kf + ifdr
}

combos <- expand.grid(ic = seq_along(c_grid),
                      ib = seq_along(beta_grid),
                      ip = seq_along(prop.diff_grid),
                      i_nsam = seq_along(nsam_grid),
                      ifdr = seq_along(fdr_grid),
                      KEEP.OUT.ATTRS = FALSE)

combos$k     <- with(combos, idx5(ic, ib, ip, i_nsam, ifdr))
combos$label <- sprintf("c=%g|beta=%g|prop.diff=%.2f|nsam=%d_%d|fdr=%.2f",
                        c_grid[combos$ic], 
                        beta_grid[combos$ib],
                        prop.diff_grid[combos$ip],
                        sapply(nsam_grid[combos$i_nsam], `[`, 1),  # grp1
                        sapply(nsam_grid[combos$i_nsam], `[`, 2),  # grp2
                        fdr_grid[combos$ifdr])


## method
methods <- c("DUET-knockoff")

## storage for each combo
fdr_store <- vector("list", K)
power_store <- vector("list", K)

for (k in seq_len(K)) {
  fdr_store[[k]] <- matrix(NA_real_, nrow = B, 
                           ncol = length(methods),
                           dimnames = list(NULL, methods))
  power_store[[k]] <- matrix(NA_real_, nrow = B, 
                             ncol = length(methods),
                             dimnames = list(NULL, methods))
}


## outer cores
n_cores_outer <- min(72, parallel::detectCores() - 1L)

## inner cores
# n_cores <- 1

## repeat seed
base_seeds <- seq(11000, 11000 + (K - 1) * 1000, by = 1000)
all_seeds  <- lapply(seq_len(K), function(k) base_seeds[k] + seq_len(B))

## ========================================================
## Flattened task parallelism (each combo × repeat is an independent task)
## ========================================================

tasks <- combos[rep(seq_len(nrow(combos)), each = B), ]
tasks$r    <- rep(seq_len(B), times = nrow(combos))
tasks$seed <- unlist(all_seeds)
stopifnot(length(unique(tasks$seed)) == nrow(tasks))

## Preload package
suppressPackageStartupMessages({
  try(library(DUETknockoff), silent = TRUE)
  try(library(dirmult), silent = TRUE)
})


# signle work function
work_fun <- function(i){
  
  set.seed(tasks$seed[i])
  
  Sys.setenv(
    MC_CORES = "1",
    OMP_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  
  ic   <- tasks$ic[i]
  ib   <- tasks$ib[i]
  ip   <- tasks$ip[i]
  i_nsam <- tasks$i_nsam[i]
  ifdr <- tasks$ifdr[i]
  k    <- tasks$k[i]
  r    <- tasks$r[i]
  
  c_val      <- c_grid[ic]
  beta_val   <- beta_grid[ib]
  prop.diff_val <- prop.diff_grid[ip]
  nsam_val   <- nsam_grid[[i_nsam]]
  fdr_target <- fdr_grid[ifdr]
  
  vals <- try(
    simulate_one_param(
      c_val      = c_val,
      beta_val   = beta_val,
      prop.diff_val = prop.diff_val,
      nsam_val   = nsam_val,
      fdr_target = fdr_target,
      seed       = tasks$seed[i]
    ),
    silent = TRUE
  )
  
  if (inherits(vals, "try-error")) {
    cat("\n--- ERROR IN SIMULATION (k =", k, ", r =", r, ") ---\n")
    print(vals)
    return(list(
      k     = k,
      r     = r,
      fdr   = rep(NA_real_, length(methods)),
      power = rep(NA_real_, length(methods))
    ))
  }
  
  return(list(
    k     = k,
    r     = r,
    fdr   = as.numeric(vals$fdr[methods]),
    power = as.numeric(vals$power[methods])
  ))
  
}

## pbmcapply
use_pb <- suppressWarnings(requireNamespace("pbmcapply", quietly = TRUE))
res_task_list <- if (use_pb) {
  pbmcapply::pbmclapply(
    X = seq_len(nrow(tasks)),
    FUN = work_fun,
    mc.cores = n_cores_outer,
    mc.preschedule = FALSE,   # 
    mc.set.seed = FALSE,
    ignore.interactive = TRUE
  )
} else {
  message("Install 'pbmcapply' to enable a progress bar: install.packages('pbmcapply')")
  mclapply(
    X = seq_len(nrow(tasks)),
    FUN = work_fun,
    mc.cores = n_cores_outer,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE,      # already manage seeds manually
    ignore.interactive = TRUE
  )
}

# Populate task-level results back into fdr_store / power_store
for (res in res_task_list) {
  fdr_store[[res$k]][res$r, ]   <- res$fdr
  power_store[[res$k]][res$r, ] <- res$power
}

## Summary output
summary_rows <- vector("list", K * length(methods))
row_i <- 0L

for (k in seq_len(K)) {
  for (m in methods) {
    row_i <- row_i + 1L
    summary_rows[[row_i]] <- data.frame(
      method     = m,
      label      = combos$label[k],
      c          = c_grid[combos$ic[k]],
      beta       = beta_grid[combos$ib[k]],
      prop.diff  = prop.diff_grid[combos$ip[k]],
      n_grp1 = nsam_grid[[combos$i_nsam[k]]][1],
      n_grp2 = nsam_grid[[combos$i_nsam[k]]][2],
      fdr_target = fdr_grid[combos$ifdr[k]],
      avg_fdr    = mean(fdr_store[[k]][, m], na.rm = TRUE),
      avg_power  = mean(power_store[[k]][, m], na.rm = TRUE),
      row.names  = NULL
    )
  }
}

res_df <- do.call(rbind, summary_rows)

print(res_df)

