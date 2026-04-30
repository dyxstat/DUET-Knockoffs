# monte carlo - CRN version with method-level seed separation
# Modified: causal OTUs, confounding, bias fixed outside Monte Carlo loop
# (following published Simulation_code_complete_overlap.R structure)
# Sensitivity version: depth_fold added as a Monte Carlo grid parameter

library(parallel)
library(MASS)
library(gtools)
library(dirmult)
library(DUETKnockoffs)
library(LOCOM)
library(readxl)
library(openxlsx)

Sys.setenv(
  MC_CORES = "1",
  OMP_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

RNGkind("Mersenne-Twister")
set.seed(20160106)

B <- 100  # number of Monte Carlo replicates

n_confounders <- 1

filter.thresh <- 0.2
n_rej_stop <- 100

bias_sd1 <- 0.5
bias_sd2 <- 0.5

depth.mu1 <- 100000

have_bias <- 1
have_diff_bias <- 1
depth.sd1 <- depth.mu1/3
depth.lower <- 10000

### read-in data (published code's template)
df <- read_xlsx('template data/throat.otu.table.xlsx')
df <- as.data.frame(df)
pi_est <- read.table('template data/fit_dirmult_pi.txt', header=FALSE, as.is=TRUE)[,1]
n.otus <- length(pi_est)

##########################################################
## Grid parameters
##########################################################

c_grid <- c(1e4, 1e5, 1e6)
beta_grid <- c(4, 6, 8, 10, 12)
prop.diff_grid <- c(0.05)
nsam_grid <- list(c(50, 50))
fdr_grid <- c(0.2)
disp_grid <- c(0.05)
depth_fold_grid <- c(50, 200, 1000)

##########################################################
## Fixed components — outside Monte Carlo loop
## (causal OTUs, confounders, bias: fixed once across all replicates)
##########################################################

## Taxa names (for result summarization)
taxa_names <- paste0("taxon", 1:n.otus)

## 1) Causal OTUs — fixed once per prop.diff
causal.otus.idx.list    <- vector("list", length(prop.diff_grid))
noncausal.otus.idx.list <- vector("list", length(prop.diff_grid))
causal.otus.list        <- vector("list", length(prop.diff_grid))
noncausal.otus.list     <- vector("list", length(prop.diff_grid))

for (ip in seq_along(prop.diff_grid)) {
  n_DA <- ceiling(n.otus * prop.diff_grid[ip])
  causal.otus.idx.list[[ip]]    <- sample(which(pi_est > 1e-8), n_DA)
  noncausal.otus.idx.list[[ip]] <- setdiff(1:n.otus, causal.otus.idx.list[[ip]])
  causal.otus.list[[ip]]        <- taxa_names[causal.otus.idx.list[[ip]]]
  noncausal.otus.list[[ip]]     <- setdiff(taxa_names, causal.otus.list[[ip]])
}

## 2) Confounding OTUs — fixed once
confounding.otus <- NULL
betaC <- 1.5
if (n_confounders > 0) {
  confounding.otus <- matrix(NA, nrow=n_confounders, ncol=5)
  w <- which(pi_est >= 0.005)
  for (r in 1:n_confounders) confounding.otus[r,] <- sort(sample(w, 5))
}

## 3) Bias factors — fixed once per prop.diff (dedicated seed, following published code)
bias.factor1.list <- vector("list", length(prop.diff_grid))
bias.factor2.list <- vector("list", length(prop.diff_grid))

for (ip in seq_along(prop.diff_grid)) {
  set.seed(0 + ip)
  bf1 <- rep(1, n.otus)
  bf2 <- rep(1, n.otus)
  if (have_bias) {
    bf1.log <- rnorm(n.otus, 0, bias_sd1)
    bf2.log <- rnorm(n.otus, 0, bias_sd2)

    if (have_diff_bias) {
      sub1 <- 1:5
      sub2 <- 11:15

      top5 <- order(pi_est, decreasing = TRUE)[1:5]
      sub1 <- c(causal.otus.idx.list[[ip]][sub1],
                setdiff(sample(noncausal.otus.idx.list[[ip]], length(noncausal.otus.idx.list[[ip]])/5), top5))
      sub2 <- c(causal.otus.idx.list[[ip]][sub2],
                setdiff(sample(noncausal.otus.idx.list[[ip]], length(noncausal.otus.idx.list[[ip]])/5), top5))

      bf1.log[sub1] <- -5
      bf2.log[sub2] <- -5
    }

    bf1 <- exp(bf1.log)
    bf2 <- exp(bf2.log)
    top_taxa <- order(pi_est, decreasing = TRUE)[1]
    bf1[top_taxa] <- 1
    bf2[top_taxa] <- 1
  }
  bias.factor1.list[[ip]] <- bf1
  bias.factor2.list[[ip]] <- bf2
}

##########################################################
## simulate_one_param: per-replicate simulation
## (causal, confounding, bias accessed from global scope — fixed)
##########################################################

simulate_one_param <- function(c_val, beta_val, nsam_val, fdr_target, disp_val, depth_fold, ip, seed) {

  ## ========== CRN: per-replicate randomness ==========
  set.seed(seed)

  n_cores <- 1

  ## Depth parameters depend on depth_fold (scenario-level)
  depth.mu2 <- depth.mu1 * depth_fold
  depth.sd2 <- depth.mu2 / 3

  ## Select causal OTUs and bias factors for this prop.diff
  causal.otus.idx <- causal.otus.idx.list[[ip]]
  causal.otus     <- causal.otus.list[[ip]]
  noncausal.otus  <- noncausal.otus.list[[ip]]
  bias.factor1    <- bias.factor1.list[[ip]]
  bias.factor2    <- bias.factor2.list[[ip]]

  n.sam.grp1 <- nsam_val[1]
  n.sam.grp2 <- nsam_val[2]
  n_sam <- n.sam.grp1 + n.sam.grp2

  ## beta.otu (deterministic given beta_val)
  beta.otu <- rep(beta_val, length(causal.otus.idx))
  beta.otu.log <- log(beta.otu)

  ## 1) Covariates C (random per replicate)
  if (n_confounders > 0) {
    C <- matrix(NA, nrow = n_sam, ncol = n_confounders)
    for (r in 1:n_confounders) {
      C[, r] <- c(runif(n.sam.grp1, -1, 1), runif(n.sam.grp2, 0, 2))
    }
  } else {
    C <- NULL
  }

  ## 2) Depth (random per replicate, rnorm)
  depth1.sim <- rnorm(n_sam, depth.mu1, depth.sd1)
  depth1.sim[depth1.sim < depth.lower] <- depth.lower
  depth1.sim <- round(depth1.sim)

  depth2.sim <- rnorm(n_sam, depth.mu2, depth.sd2)
  depth2.sim[depth2.sim < depth.lower] <- depth.lower
  depth2.sim <- round(depth2.sim)

  ## 3) Group labels
  Y <- c(rep(0, n.sam.grp1), rep(1, n.sam.grp2))

  ## 4) Dirichlet sampling (sample-level variation)
  alpha <- c_val * pi_est
  pi.table.sim <- rdirichlet(n = n_sam, alpha)
  rownames(pi.table.sim) <- paste0("sub", 1:n_sam)
  colnames(pi.table.sim) <- taxa_names

  pi.table1.sim <- pi.table.sim
  pi.table2.sim <- pi.table.sim

  ## 5) Causal effect (deterministic, using fixed causal.otus.idx)
  pi.table1.sim[, causal.otus.idx] <- pi.table1.sim[, causal.otus.idx] * exp(Y %*% t(beta.otu.log))
  pi.table2.sim[, causal.otus.idx] <- pi.table2.sim[, causal.otus.idx] * exp(Y %*% t(beta.otu.log))

  ## 6) Confounding effect (deterministic, using fixed confounding.otus)
  if (n_confounders > 0) {
    for (r in 1:n_confounders) {
      conf_factor <- betaC^C[, r]
      pi.table1.sim[, confounding.otus[r,]] <- sweep(pi.table1.sim[, confounding.otus[r,], drop=FALSE], 1, conf_factor, "*")
      pi.table2.sim[, confounding.otus[r,]] <- sweep(pi.table2.sim[, confounding.otus[r,], drop=FALSE], 1, conf_factor, "*")
    }
  }

  ## 7) Bias (deterministic, using fixed bias.factor)
  if (have_bias == 1) {
    pi.table1.sim <- sweep(pi.table1.sim, 2, bias.factor1, "*")
    pi.table2.sim <- sweep(pi.table2.sim, 2, bias.factor2, "*")
  }

  ## 8) Normalizing
  pi.table1.sim <- pi.table1.sim / rowSums(pi.table1.sim)
  pi.table2.sim <- pi.table2.sim / rowSums(pi.table2.sim)

  ## 9) Count data (Dirichlet-Multinomial or Multinomial)
  otu.table1.sim <- matrix(0, n_sam, n.otus)
  otu.table2.sim <- matrix(0, n_sam, n.otus)
  colnames(otu.table1.sim) <- colnames(otu.table2.sim) <- taxa_names

  for (i in 1:n_sam) {
    if (disp_val < 1e-8) {
      otu.table1.sim[i,] <- rmultinom(1, depth1.sim[i], pi.table1.sim[i,])
      otu.table2.sim[i,] <- rmultinom(1, depth2.sim[i], pi.table2.sim[i,])
    } else {
      otu.table1.sim[i,] <- simPop(J = 1, n = depth1.sim[i], pi = pi.table1.sim[i,], theta = disp_val)$data
      otu.table2.sim[i,] <- simPop(J = 1, n = depth2.sim[i], pi = pi.table2.sim[i,], theta = disp_val)$data
    }
  }

  rn <- sprintf("id_%03d", seq_len(nrow(otu.table1.sim)))
  rownames(otu.table1.sim) <- rn
  rownames(otu.table2.sim) <- rn

  ## 10) Filter & pool
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

  ## Safety check
  if (length(otus.keep.pool) == 0 || ncol(otu.table.sim.pool.filter) == 0) {
    fdr_vec   <- c(`DUET-Knockoffs`=NA, Com2seq=NA, Locom_Com_P=NA, Locom_Com_Count=NA,
                   Locom_16s=NA, Locom_shotgun=NA)
    power_vec <- fdr_vec
    return(list(fdr = fdr_vec, power = power_vec))
  }

  ## ==== Prepare method inputs ====
  W1 <- otu.table1.sim[, otus.keep.pool, drop = FALSE]
  W2 <- otu.table2.sim[, otus.keep.pool, drop = FALSE]
  M1 <- rowSums(W1)
  M2 <- rowSums(W2)

  # data_x_1 <- data_x_2 <- as.data.frame(cbind(C[,1], C[,2], C[,3]))
  # colnames(data_x_1) <- c("X_1","X_2","X_3")
  # colnames(data_x_2) <- c("X_1","X_2","X_3")

  data_x_1 <- data_x_2 <- as.data.frame(C[,1])
  colnames(data_x_1) <- "X_1"
  colnames(data_x_2) <- "X_1"

  count_16s <- cbind(data_x_1, M = M1, W1)
  count_SMS <- cbind(data_x_2, M = M2, W2)
  count_SK <- rbind(count_16s, count_SMS)

  class_K <- factor(rep(1:2, each = n_sam))
  y <- rep(c(rep(1, n.sam.grp1), rep(2, n.sam.grp2)), times = 2)

  n_x1 <- 1
  # data_x <- as.data.frame(count_SK[, 1:n_x1])
  data_x <- as.data.frame(count_SK[, 1:n_x1, drop = FALSE])
  W <- as.data.frame(count_SK[, -c(1:(n_x1 + 1))])
  M <- count_SK[, n_x1 + 1]

  T_var <- which(colnames(W) %in% causal.otus)

  ## ========== Method calls (beta-dependent seed) ==========

  beta_int <- as.integer(beta_val)

  ## DUETKnockoffs
  res.DUETKnockoffs <- DUET_Knockoffs(
    W = W, class_K = class_K, data_x = data_x, M = M, y = y, T_var = T_var,
    fdr = fdr_target, test_statistic = "DE", filter_statistics = 3, test1 = "wilcox.test",
    offset = 1
  )
  duetknockoffs.fdr <- unname(res.DUETKnockoffs$res[2])
  duetknockoffs.power <- unname(res.DUETKnockoffs$res[3])

  ## Com2seq
  Y1 <- matrix(as.integer(Y==1), ncol=1)
  Y2 <- Y1
  C1 <- as.matrix(C)
  C2 <- as.matrix(C)

  res.Com2seq <- Com2seq(
    table1 = W1, table2 = W2,
    Y1 = Y1, Y2 = Y2,
    C1 = C1, C2 = C2,
    fdr.nominal = fdr_target, n.cores = n_cores,
    n.perm.max = 50000, n.rej.stop = n_rej_stop,
    seed = seed + 500000L + beta_int * 1000L
  )

  ## LOCOM (using beta-dependent seed)
  res.locom1 <- locom(otu.table = otu.table1.sim.filter, Y = Y, C = C,
                      n.perm.max = 50000, fdr.nominal = fdr_target,
                      n.cores = n_cores, n.rej.stop = n_rej_stop)

  res.locom2 <- locom(otu.table = otu.table2.sim.filter, Y = Y, C = C,
                      n.perm.max = 50000, fdr.nominal = fdr_target,
                      n.cores = n_cores, n.rej.stop = n_rej_stop)

  res.locom.pool <- locom(otu.table = otu.table.sim.pool.filter, Y = Y, C = C,
                          n.perm.max = 50000, fdr.nominal = fdr_target,
                          n.cores = n_cores, n.rej.stop = n_rej_stop)

  ## p.combine (Locom 16s + shotgun)
  name1 <- colnames(res.locom1$p.otu)
  name2 <- colnames(res.locom2$p.otu)
  common.name <- intersect(name1, name2)
  j.mat1 <- match(common.name, name1)
  j.mat2 <- match(common.name, name2)

  p.both <- rbind(res.locom1$p.otu[j.mat1], res.locom2$p.otu[j.mat2])
  p.comp <- pcauchy(apply(tan((0.5 - p.both)*pi), 2, mean), lower.tail = FALSE)

  p.comp.name <- common.name
  if (length(res.locom1$p.otu[-j.mat1]) > 0) {
    p.comp <- c(p.comp, res.locom1$p.otu[-j.mat1])
    p.comp.name <- c(p.comp.name, name1[-j.mat1])
  }
  if (length(res.locom2$p.otu[-j.mat2]) > 0) {
    p.comp <- c(p.comp, res.locom2$p.otu[-j.mat2])
    p.comp.name <- c(p.comp.name, name2[-j.mat2])
  }

  p.comp <- matrix(p.comp, nrow=1)
  q.comp <- matrix(p.adjust(p.comp, method="BH"), nrow=1)
  colnames(q.comp) <- p.comp.name

  ## summarize results
  summarize_otu_results <- function(qvalue, causal.otus, noncausal.otus, fdr.target=fdr_target) {
    otu.detected <- colnames(qvalue)[which(qvalue < fdr.target)]
    n.otu <- length(otu.detected)
    if (n.otu > 0) {
      sen <- sum(otu.detected %in% causal.otus) / length(causal.otus)
      fdr <- (n.otu - sum(otu.detected %in% causal.otus)) / n.otu
    } else {
      sen <- 0
      fdr <- 0
    }
    list(n.otu=n.otu, sen=sen, fdr=fdr)
  }

  otu.new.omni <- summarize_otu_results(res.Com2seq$q.taxa.omni, causal.otus, noncausal.otus)
  otu.comp.locom <- summarize_otu_results(q.comp, causal.otus, noncausal.otus)
  otu.pool.locom <- summarize_otu_results(res.locom.pool$q.otu, causal.otus, noncausal.otus)
  otu.locom.1 <- summarize_otu_results(res.locom1$q.otu, causal.otus, noncausal.otus)
  otu.locom.2 <- summarize_otu_results(res.locom2$q.otu, causal.otus, noncausal.otus)

  fdr_vec <- c(`DUET-Knockoffs`=duetknockoffs.fdr, Com2seq=otu.new.omni$fdr,
               Locom_Com_P=otu.comp.locom$fdr, Locom_Com_Count=otu.pool.locom$fdr,
               Locom_16s=otu.locom.1$fdr, Locom_shotgun=otu.locom.2$fdr)

  power_vec <- c(`DUET-Knockoffs`=duetknockoffs.power, Com2seq=otu.new.omni$sen,
                 Locom_Com_P=otu.comp.locom$sen, Locom_Com_Count=otu.pool.locom$sen,
                 Locom_16s=otu.locom.1$sen, Locom_shotgun=otu.locom.2$sen)

  list(fdr = fdr_vec, power = power_vec)
}


## -----------------------------
## CRN Monte Carlo over grid
## -----------------------------

methods <- c("DUET-Knockoffs", "Com2seq", "Locom_Com_P", "Locom_Com_Count", "Locom_16s", "Locom_shotgun")

Kc  <- length(c_grid)
Kp  <- length(prop.diff_grid)
Kn  <- length(nsam_grid)
Kf  <- length(fdr_grid)
Kb  <- length(beta_grid)
Kd  <- length(disp_grid)
Kdf <- length(depth_fold_grid)

## combos without beta (but includes disp and depth_fold)
combos0 <- expand.grid(
  ic = seq_along(c_grid),
  ip = seq_along(prop.diff_grid),
  i_nsam = seq_along(nsam_grid),
  ifdr = seq_along(fdr_grid),
  idisp = seq_along(disp_grid),
  idf = seq_along(depth_fold_grid),
  KEEP.OUT.ATTRS = FALSE
)
K0 <- nrow(combos0)
combos0$k0 <- seq_len(K0)
combos0$label <- sprintf("c=%g|prop.diff=%.2f|nsam=%d_%d|fdr=%.2f|disp=%.3f|depth_fold=%g",
                         c_grid[combos0$ic],
                         prop.diff_grid[combos0$ip],
                         sapply(nsam_grid[combos0$i_nsam], `[`, 1),
                         sapply(nsam_grid[combos0$i_nsam], `[`, 2),
                         fdr_grid[combos0$ifdr],
                         disp_grid[combos0$idisp],
                         depth_fold_grid[combos0$idf])

## store: each combo -> B x beta x methods
fdr_store <- vector("list", K0)
power_store <- vector("list", K0)
for (k0 in seq_len(K0)) {
  fdr_store[[k0]] <- array(NA_real_, dim = c(B, Kb, length(methods)),
                           dimnames = list(NULL, as.character(beta_grid), methods))
  power_store[[k0]] <- array(NA_real_, dim = c(B, Kb, length(methods)),
                             dimnames = list(NULL, as.character(beta_grid), methods))
}

## tasks: combos0 x B
tasks <- combos0[rep(seq_len(K0), each = B), ]
tasks$r <- rep(seq_len(B), times = K0)

## seed stream depends ONLY on (k0, r) — not beta
## k0 already incorporates idf (via expand.grid), so distinct depth_fold scenarios get distinct seeds
base_seeds <- seq(1000000, by = 3000, length.out = K0)
tasks$seed <- base_seeds[rep(seq_len(K0), each = B)] + tasks$r
stopifnot(length(unique(tasks$seed)) == nrow(tasks))

## cores
n_cores_outer <- min(180, parallel::detectCores() - 1L)

suppressPackageStartupMessages({
  try(library(DUETKnockoffs), silent = TRUE)
  try(library(LOCOM), silent = TRUE)
  try(library(dirmult), silent = TRUE)
})

## worker
work_fun_crn <- function(i){

  Sys.setenv(
    MC_CORES="1",
    OMP_NUM_THREADS="1",
    MKL_NUM_THREADS="1",
    OPENBLAS_NUM_THREADS="1",
    NUMEXPR_NUM_THREADS="1"
  )

  ic <- tasks$ic[i]
  ip <- tasks$ip[i]
  i_nsam <- tasks$i_nsam[i]
  ifdr <- tasks$ifdr[i]
  idisp <- tasks$idisp[i]
  idf <- tasks$idf[i]
  k0 <- tasks$k0[i]
  r <- tasks$r[i]
  seed_r <- tasks$seed[i]

  c_val <- c_grid[ic]
  nsam_val <- nsam_grid[[i_nsam]]
  fdr_target <- fdr_grid[ifdr]
  disp_val <- disp_grid[idisp]
  depth_fold_val <- depth_fold_grid[idf]

  ## For THIS replicate seed_r, run ALL betas (CRN)
  out_beta <- lapply(seq_along(beta_grid), function(ib){
    beta_val <- beta_grid[ib]
    vals <- try(
      simulate_one_param(
        c_val=c_val, beta_val=beta_val,
        nsam_val=nsam_val,
        fdr_target=fdr_target, disp_val=disp_val,
        depth_fold=depth_fold_val,
        ip=ip, seed=seed_r
      ),
      silent = TRUE
    )
    if (inherits(vals, "try-error")) {
      return(list(
        fdr = setNames(rep(NA_real_, length(methods)), methods),
        power = setNames(rep(NA_real_, length(methods)), methods)
      ))
    }
    list(fdr = vals$fdr, power = vals$power)
  })

  fdr_mat <- do.call(rbind, lapply(out_beta, function(z) as.numeric(z$fdr[methods])))
  pow_mat <- do.call(rbind, lapply(out_beta, function(z) as.numeric(z$power[methods])))

  rownames(fdr_mat) <- as.character(beta_grid)
  rownames(pow_mat) <- as.character(beta_grid)
  colnames(fdr_mat) <- methods
  colnames(pow_mat) <- methods

  list(k0=k0, r=r, fdr_mat=fdr_mat, pow_mat=pow_mat)
}

use_pb <- suppressWarnings(requireNamespace("pbmcapply", quietly = TRUE))
res_task_list <- if (use_pb) {
  pbmcapply::pbmclapply(
    X = seq_len(nrow(tasks)),
    FUN = work_fun_crn,
    mc.cores = n_cores_outer,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE,
    ignore.interactive = TRUE
  )
} else {
  mclapply(
    X = seq_len(nrow(tasks)),
    FUN = work_fun_crn,
    mc.cores = n_cores_outer,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE,
    ignore.interactive = TRUE
  )
}

## fill back
for (res in res_task_list) {
  if (is.null(res) || is.null(res$k0)) next
  fdr_store[[res$k0]][res$r, , ] <- res$fdr_mat
  power_store[[res$k0]][res$r, , ] <- res$pow_mat
}

## summarize to long df
summary_rows <- vector("list", K0 * Kb * length(methods))
row_i <- 0L

for (k0 in seq_len(K0)) {
  for (ib in seq_along(beta_grid)) {
    for (m in methods) {
      row_i <- row_i + 1L
      summary_rows[[row_i]] <- data.frame(
        method = m,
        label  = combos0$label[k0],
        c      = c_grid[combos0$ic[k0]],
        beta   = beta_grid[ib],
        prop.diff = prop.diff_grid[combos0$ip[k0]],
        n_grp1 = nsam_grid[[combos0$i_nsam[k0]]][1],
        n_grp2 = nsam_grid[[combos0$i_nsam[k0]]][2],
        fdr_target = fdr_grid[combos0$ifdr[k0]],
        disp   = disp_grid[combos0$idisp[k0]],
        depth_fold = depth_fold_grid[combos0$idf[k0]],
        avg_fdr   = mean(fdr_store[[k0]][, ib, m], na.rm = TRUE),
        avg_power = mean(power_store[[k0]][, ib, m], na.rm = TRUE),
        row.names = NULL
      )
    }
  }
}

res_df <- do.call(rbind, summary_rows)
print(res_df)

#############################################
## END
#############################################
