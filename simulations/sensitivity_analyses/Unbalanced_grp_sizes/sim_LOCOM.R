# monte carlo 100 - CRN version with method-level seed separation
# LOCOM methods version

library(parallel)
library(MASS)
library(gtools)
library(dirmult)
library(LOCOM)
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

df <- readRDS("ZhuF_count.rds") # YeZ_count, Dha_count
df <- df[, colSums(df) > 0]

pi_est <- DM.MoM(df)$pi
n.otus <- length(pi_est)

## -----------------------------
## 100 times Monte Carlo over grid (LOCOM_16s, LOCOM_shotgun, LOCOM_Com_p, LOCOM_Com_count)
## -----------------------------

c_grid <- c(1e3, 1e4, 1e5, 1e6)
beta_grid <- c(5, 6, 7, 8, 9, 10, 11, 12)
prop.diff_grid <- c(0.15)
nsam_grid <- list(c(20, 30), c(40, 60))
fdr_grid <- c(0.2)

methods <- c("Locom_Com_P", "Locom_Com_Count", "Locom_16s", "Locom_shotgun")

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
## single simulation
## -----------------------------
simulate_one_param <- function(c_val, beta_val, prop.diff_val, nsam_val, fdr_target, seed, ib) {

  set.seed(seed)

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
  bias.factor1 <- rep(1, n.otus)
  bias.factor2 <- rep(1, n.otus)
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
  c <- c_val    # 1e4, 1e5, 1e6
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

  pi.table1.sim <- pi.table1.sim / rowSums(pi.table1.sim)
  pi.table2.sim <- pi.table2.sim / rowSums(pi.table2.sim)

  ################################################################

  # 7) Poisson-Gamma sampling
  mu1 <- pi.table1.sim * depth1.sim    # mean absolute abundance matrix for 16s
  mu2 <- pi.table2.sim * depth2.sim

  phi_safe <- pmax(phi, 1e-8)

  otu.table1.sim <- matrix(0, n_sam, n.otus)
  otu.table2.sim <- matrix(0, n_sam, n.otus)

  idx <- (phi >= 1e-8)     # Poisson-gamma for phi != 0
  idx_p  <- !idx           # Poisson sampling for phi = 0

  # 16S count matrix
  if (any(idx)) {
    mu1_pg <- mu1[, idx, drop = FALSE]
    A1 <- rgamma(n = n_sam * sum(idx),
                 shape = rep(1 / phi_safe[idx], each = n_sam),
                 rate = rep(1 / phi_safe[idx], each = n_sam) / c(mu1_pg))
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
    mu2_pg <- mu2[, idx, drop = FALSE]
    A2 <- rgamma(n = n_sam * sum(idx),
                 shape = rep(1 / phi_safe[idx], each = n_sam),
                 rate = rep(1 / phi_safe[idx], each = n_sam) / c(mu2_pg))
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

  # ==== LOCOM ==== #
  set.seed(seed + 100000L + ib * 1000L)
  res.locom1 <- locom(otu.table = otu.table1.sim.filter, Y = Y, C = C,
                      n.perm.max = 50000, fdr.nominal = fdr_target,
                      n.cores = n_cores, n.rej.stop = n_rej_stop)

  set.seed(seed + 200000L + ib * 1000L)
  res.locom2 <- locom(otu.table = otu.table2.sim.filter, Y = Y, C = C,
                      n.perm.max = 50000, fdr.nominal = fdr_target,
                      n.cores = n_cores, n.rej.stop = n_rej_stop)

  set.seed(seed + 300000L + ib * 1000L)
  res.locom.pool <- locom(otu.table = otu.table.sim.pool.filter, Y = Y, C = C,
                          n.perm.max = 50000, fdr.nominal = fdr_target,
                          n.cores = n_cores, n.rej.stop = n_rej_stop)

  ###########
  # p.combine
  ###########

  name1 <- colnames(res.locom1$p.otu)    # res.locom1 otu name
  name2 <- colnames(res.locom2$p.otu)    # res.locom2 otu name
  common.name <- intersect(name1, name2) # common otu name
  j.mat1 <- match(common.name, name1)  # index for otu name in name1
  j.mat2 <- match(common.name, name2)  # index for otu name in name2

  p.both <- rbind(res.locom1$p.otu[j.mat1],
                  res.locom2$p.otu[j.mat2])  # rbind common p-value
  p.comp <- pcauchy(apply( tan( (0.5 - p.both)*pi ), 2, mean),
                    lower.tail = F)  # use pcauchy to combine p-value

  p.comp.name <- common.name  #

  # combine single p-values from res.locom1 and res.locom2 into p.comp
  if (length(res.locom1$p.otu[-j.mat1]) > 0) {
    p.comp <- c(p.comp, res.locom1$p.otu[-j.mat1])
    p.comp.name <- c(p.comp.name, name1[-j.mat1])
  }
  if (length(res.locom2$p.otu[-j.mat2]) > 0) {
    p.comp <- c(p.comp, res.locom2$p.otu[-j.mat2])
    p.comp.name <- c(p.comp.name, name2[-j.mat2])
  }  #

  # Benjamini-Hochberg adjust for p.comp
  p.comp <- matrix(p.comp, nrow=1)
  q.comp<- matrix(p.adjust(p.comp, method ="BH"), nrow=1)
  colnames(p.comp) <- p.comp.name
  colnames(q.comp) <- p.comp.name

  # -----------------------------------
  # Summarizing results
  # -----------------------------------

  ## define function, out:n.out, se, sep, fdr
  summarize_otu_results <- function(qvalue, causal.otus, noncausal.otus, fdr.target=fdr_target) {

    otu.detected = colnames(qvalue)[which(qvalue < fdr.target)]
    n.otu = length(otu.detected)

    if (n.otu > 0) {
      sen = sum(otu.detected %in% causal.otus)/length(causal.otus)
      fdr = n.otu - sum(otu.detected %in% causal.otus)
      fdr = fdr/n.otu
    } else {
      sen = 0
      fdr = 0
    }

    out = list(n.otu=n.otu, sen=sen, fdr=fdr)

    return(out)
  }

  # Benchmark q-value summary
  otu.comp.locom <- summarize_otu_results(q.comp, causal.otus, noncausal.otus)
  otu.pool.locom <- summarize_otu_results(res.locom.pool$q.otu, causal.otus, noncausal.otus)
  otu.locom.1 <- summarize_otu_results(res.locom1$q.otu, causal.otus, noncausal.otus)
  otu.locom.2 <- summarize_otu_results(res.locom2$q.otu, causal.otus, noncausal.otus)

  fdr_vec <- c(Locom_Com_P=otu.comp.locom$fdr, Locom_Com_Count=otu.pool.locom$fdr,
               Locom_16s=otu.locom.1$fdr, Locom_shotgun=otu.locom.2$fdr)

  power_vec <- c(Locom_Com_P=otu.comp.locom$sen, Locom_Com_Count=otu.pool.locom$sen,
                 Locom_16s=otu.locom.1$sen, Locom_shotgun=otu.locom.2$sen)

  list(fdr = fdr_vec, power = power_vec)
}


## -----------------------------
## CRN Monte Carlo over grid (beta NOT in combos)
## -----------------------------

Kc <- length(c_grid)
Kp <- length(prop.diff_grid)
Kn <- length(nsam_grid)
Kf <- length(fdr_grid)
Kb <- length(beta_grid)

## combos without beta
combos0 <- expand.grid(
  ic = seq_along(c_grid),
  ip = seq_along(prop.diff_grid),
  i_nsam = seq_along(nsam_grid),
  ifdr = seq_along(fdr_grid),
  KEEP.OUT.ATTRS = FALSE
)
K0 <- nrow(combos0)
combos0$k0 <- seq_len(K0)
combos0$label <- sprintf("c=%g|prop.diff=%.2f|nsam=%d_%d|fdr=%.2f",
                         c_grid[combos0$ic],
                         prop.diff_grid[combos0$ip],
                         sapply(nsam_grid[combos0$i_nsam], `[`, 1),
                         sapply(nsam_grid[combos0$i_nsam], `[`, 2),
                         fdr_grid[combos0$ifdr])

## store: each combo -> B × beta × methods
fdr_store <- vector("list", K0)
power_store <- vector("list", K0)
for (k0 in seq_len(K0)) {
  fdr_store[[k0]] <- array(NA_real_, dim = c(B, Kb, length(methods)),
                           dimnames = list(NULL, as.character(beta_grid), methods))
  power_store[[k0]] <- array(NA_real_, dim = c(B, Kb, length(methods)),
                             dimnames = list(NULL, as.character(beta_grid), methods))
}

## tasks: combos0 × B
tasks <- combos0[rep(seq_len(K0), each = B), ]
tasks$r <- rep(seq_len(B), times = K0)

## seed stream depends ONLY on (k0, r) — not beta
base_seeds <- seq(1000000, by = 5000, length.out = K0)
tasks$seed <- base_seeds[rep(seq_len(K0), each = B)] + tasks$r
stopifnot(length(unique(tasks$seed)) == nrow(tasks))

## cores
n_cores_outer <- min(72, parallel::detectCores() - 1L)

## Preload package
suppressPackageStartupMessages({
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
  k0 <- tasks$k0[i]
  r <- tasks$r[i]
  seed_r <- tasks$seed[i]

  c_val <- c_grid[ic]
  prop.diff_val <- prop.diff_grid[ip]
  nsam_val <- nsam_grid[[i_nsam]]
  fdr_target <- fdr_grid[ifdr]

  ## For THIS replicate seed_r, run ALL betas (CRN)
  out_beta <- lapply(seq_along(beta_grid), function(ib){
    beta_val <- beta_grid[ib]
    vals <- try(
      simulate_one_param(
        c_val=c_val, beta_val=beta_val,
        prop.diff_val=prop.diff_val, nsam_val=nsam_val,
        fdr_target=fdr_target, seed=seed_r, ib=ib
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

## pbmcapply
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
  message("Install 'pbmcapply' to enable a progress bar: install.packages('pbmcapply')")
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
