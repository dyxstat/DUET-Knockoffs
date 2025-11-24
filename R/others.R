#' @keywords internal
is_posdef <- function (A, tol = 1e-09)
{
  p = nrow(matrix(A))
  if (p < 500) {
    lambda_min = min(eigen(A)$values)
  }
  else {
    oldw <- getOption("warn")
    options(warn = -1)
    lambda_min = RSpectra::eigs(A, 1, which = "SM", opts = list(retvec = FALSE,
                                                                maxitr = 100, tol))$values
    options(warn = oldw)
    if (length(lambda_min) == 0) {
      lambda_min = min(eigen(A)$values)
    }
  }
  return(lambda_min > tol * 10)
}

#' @keywords internal
create.fixed <- function (X, method = c("sdp", "equi"), sigma = NULL, y = NULL,
                          randomize = F)
{
  method = match.arg(method)
  n = nrow(X)
  p = ncol(X)
  if (n <= p)
    stop("Input X must have dimensions n > p")
  else if (n < 2 * p) {
    warning("Input X has dimensions p < n < 2p. ", "Augmenting the model with extra rows.",
            immediate. = T)
    X.svd = svd(X, nu = n, nv = 0)
    u2 = X.svd$u[, (p + 1):n]
    X = rbind(X, matrix(0, 2 * p - n, p))
    if (is.null(sigma)) {
      if (is.null(y)) {
        stop("Either the noise level \"sigma\" or the response variables \"y\" must\n             be provided in order to augment the data with extra rows.")
      }
      else {
        sigma = sqrt(mean((t(u2) %*% y)^2))
      }
    }
    if (randomize)
      y.extra = rnorm(2 * p - n, sd = sigma)
    else y.extra = with_seed(0, rnorm(2 * p - n, sd = sigma))
    y = c(y, y.extra)
  }
  X = normc(X, center = F)
  Xk = switch(match.arg(method), equi = create_equicorrelated(X,
                                                              randomize), sdp = create_sdp(X, randomize))
  structure(list(X = X, Xk = Xk, y = y), class = "knockoff.variables")
}

#' @keywords internal
with_seed <- function (seed, expr)
{
  seed.old = if (exists(".Random.seed"))
    .Random.seed
  else NULL
  set.seed(seed)
  on.exit({
    if (is.null(seed.old)) {
      if (exists(".Random.seed")) rm(.Random.seed, envir = .GlobalEnv)
    } else {
      .Random.seed <<- seed.old
    }
  })
  expr
}

#' @keywords internal
normc <- function (X, center = T)
{
  X.centered = scale(X, center = center, scale = F)
  X.scaled = scale(X.centered, center = F, scale = sqrt(colSums(X.centered^2)))
  X.scaled[, ]
}

#' @keywords internal
create_equicorrelated <- function (X, randomize)
{
  X.svd = decompose(X, randomize)
  if (any(X.svd$d <= 1e-05 * max(X.svd$d)))
    stop(paste("Data matrix is rank deficient.", "Equicorrelated knockoffs will have no power."))
  lambda_min = min(X.svd$d)^2
  s = min(2 * lambda_min, 1)
  s_diff = pmax(0, 2 * s - (s/X.svd$d)^2)
  X_ko = (X.svd$u %*diag% (X.svd$d - s/X.svd$d) + X.svd$u_perp %*diag%
            sqrt(s_diff)) %*% t(X.svd$v)
}

#' @keywords internal
decompose <- function (X, randomize)
{
  n = nrow(X)
  p = ncol(X)
  stopifnot(n >= 2 * p)
  result = canonical_svd(X)
  Q = qr.Q(qr(cbind(result$u, matrix(0, n, p))))
  u_perp = Q[, (p + 1):(2 * p)]
  if (randomize) {
    Q = qr.Q(qr(rnorm_matrix(p, p)))
    u_perp = u_perp %*% Q
  }
  result$u_perp = u_perp
  result
}

#' @keywords internal
canonical_svd <- function (X)
{
  X.svd = tryCatch({
    svd(X)
  }, warning = function(w) {
  }, error = function(e) {
    stop("SVD failed in the creation of fixed-design knockoffs. Try upgrading R to version >= 3.3.0")
  }, finally = {
  })
  for (j in 1:min(dim(X))) {
    i = which.max(abs(X.svd$u[, j]))
    if (X.svd$u[i, j] < 0) {
      X.svd$u[, j] = -X.svd$u[, j]
      X.svd$v[, j] = -X.svd$v[, j]
    }
  }
  return(X.svd)
}

#' @keywords internal
rnorm_matrix <- function (n, p, mean = 0, sd = 1)
{
  matrix(rnorm(n * p, mean, sd), nrow = n, ncol = p)
}

#' @keywords internal
`%*diag%` <- function(X, d) {
  t(t(X) * d)
}

#' @keywords internal
`%diag*%` <- function(d, X) {
  d * X
}

#' @keywords internal
divide.sdp <- function (Sigma, max.size)
{
  p = ncol(Sigma)
  Eps = matrix(rnorm(p * p), p) * 1e-06
  dissimilarity = 1 - abs(cov2cor(Sigma) + Eps)
  distance = as.dist(dissimilarity)
  fit = hclust(distance, method = "single")
  n.blocks.min = 1
  n.blocks.max = ncol(Sigma)
  for (it in 1:100) {
    n.blocks = ceiling((n.blocks.min + n.blocks.max)/2)
    clusters = cutree(fit, k = n.blocks)
    size = max(table(clusters))
    if (size <= max.size) {
      n.blocks.max = n.blocks
    }
    if (size >= max.size) {
      n.blocks.min = n.blocks
    }
    if (n.blocks.min == n.blocks.max) {
      break
    }
  }
  clusters.new = merge.clusters(clusters, max.size)
  while (sum(clusters.new != clusters) > 0) {
    clusters = clusters.new
    clusters.new = merge.clusters(clusters, max.size)
  }
  clusters = clusters.new
  subSigma = vector("list", max(clusters))
  for (k in 1:length(subSigma)) {
    indices_k = clusters == k
    subSigma[[k]] = Sigma[indices_k, indices_k]
  }
  structure(list(clusters = clusters, subSigma = subSigma),
            class = "knockoff.clusteredCovariance")
}

#' @keywords internal
merge.clusters <- function (clusters, max.size)
{
  cluster.sizes = table(clusters)
  clusters.new = rep(0, length(clusters))
  g = 1
  g.size = 0
  for (k in 1:max(clusters)) {
    if (g.size + cluster.sizes[k] > max.size) {
      g = g + 1
      g.size = 0
    }
    clusters.new[clusters == k] = g
    g.size = g.size + cluster.sizes[k]
  }
  return(clusters.new)
}

#' @keywords internal
create_sdp <- function (X, randomize)
{
  X.svd = decompose(X, randomize)
  tol = 1e-05
  d = X.svd$d
  d_inv = 1/d
  d_zeros = d <= tol * max(d)
  if (any(d_zeros)) {
    warning(paste("Data matrix is rank deficient.", "Model is not identifiable, but proceeding with SDP knockoffs"),
            immediate. = T)
    d_inv[d_zeros] = 0
  }
  G = (X.svd$v %*diag% d^2) %*% t(X.svd$v)
  G_inv = (X.svd$v %*diag% d_inv^2) %*% t(X.svd$v)
  s = create.solve_sdp(G)
  s[s <= tol] = 0
  C.svd = canonical_svd(2 * diag(s) - (s %diag*% G_inv %*diag%
                                         s))
  X_ko = X - (X %*% G_inv %*diag% s) + (X.svd$u_perp %*diag%
                                          sqrt(pmax(0, C.svd$d))) %*% t(C.svd$v)
}

#'create.gaussian
#'
#' @export
create.gaussian <- function (X, mu, Sigma, method = c("asdp", "sdp", "equi"), diag_s = NULL)
{
  method = match.arg(method)
  if ((nrow(Sigma) <= 500) && method == "asdp") {
    method = "sdp"
  }
  if (is.null(diag_s)) {
    diag_s = diag(switch(match.arg(method), equi = create.solve_equi(Sigma),
                         sdp = create.solve_sdp(Sigma), asdp = create.solve_asdp(Sigma)))
  }
  if (is.null(dim(diag_s))) {
    diag_s = diag(diag_s, length(diag_s))
  }
  if (all(diag_s == 0)) {
    warning("The conditional knockoff covariance matrix is not positive definite. Knockoffs will have no power.")
    return(X)
  }
  SigmaInv_s = solve(Sigma, diag_s)
  mu_k = X - sweep(X, 2, mu, "-") %*% SigmaInv_s
  Sigma_k = 2 * diag_s - diag_s %*% SigmaInv_s
  X_k = mu_k + matrix(rnorm(ncol(X) * nrow(X)), nrow(X)) %*%
    chol(Sigma_k)
}

#' @keywords internal
create.solve_equi <- function (Sigma)
{
  stopifnot(isSymmetric(Sigma))
  p = nrow(Sigma)
  tol = 1e-10
  G = stats::cov2cor(Sigma)
  if (!is_posdef(G)) {
    stop("The covariance matrix is not positive-definite: cannot solve SDP",
         immediate. = T)
  }
  if (p > 2) {
    converged = FALSE
    maxitr = 10000
    while (!converged) {
      lambda_min = RSpectra::eigs(G, 1, which = "SR",
                                  opts = list(retvec = FALSE, maxitr = 1e+05,
                                              tol = 1e-08))$values
      if (length(lambda_min) == 1) {
        converged = TRUE
      }
      else {
        if (maxitr > 1e+08) {
          warning("In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the \n                covariance matrix. RSpectra::eigs did not converge. Giving up and computing full SVD with built-in R function.",
                  immediate. = T)
          lambda_min = eigen(G, symmetric = T, only.values = T)$values[p]
          converged = TRUE
        }
        else {
          warning("In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the \n                covariance matrix. RSpectra::eigs did not converge. Trying again with increased number of iterations.",
                  immediate. = T)
          maxitr = maxitr * 10
        }
      }
    }
  }
  else {
    lambda_min = eigen(G, symmetric = T, only.values = T)$values[p]
  }
  if (lambda_min < 0) {
    stop("In creation of equi-correlated knockoffs, while computing the smallest eigenvalue of the \n                covariance matrix. The covariance matrix is not positive-definite.")
  }
  s = rep(1, nrow(Sigma)) * min(2 * lambda_min, 1)
  psd = 0
  s_eps = 1e-08
  while (psd == 0) {
    psd = is_posdef(2 * G - diag(s * (1 - s_eps), length(s)))
    if (!psd) {
      s_eps = s_eps * 10
    }
  }
  s = s * (1 - s_eps)
  return(s * diag(Sigma))
}

#' @keywords internal
create.solve_sdp <- function (Sigma, gaptol = 1e-06, maxit = 1000, verbose = FALSE)
{
  stopifnot(isSymmetric(Sigma))
  G = cov2cor(Sigma)
  p = dim(G)[1]
  if (!is_posdef(G)) {
    warning("The covariance matrix is not positive-definite: knockoffs may not have power.",
            immediate. = T)
  }
  Cl1 = rep(0, p)
  Al1 = -Matrix::Diagonal(p)
  Cl2 = rep(1, p)
  Al2 = Matrix::Diagonal(p)
  d_As = c(diag(p))
  As = Matrix::Diagonal(length(d_As), x = d_As)
  As = As[which(Matrix::rowSums(As) > 0), ]
  Cs = c(2 * G)
  A = cbind(Al1, Al2, As)
  C = matrix(c(Cl1, Cl2, Cs), 1)
  K = NULL
  K$s = p
  K$l = 2 * p
  b = rep(1, p)
  OPTIONS = NULL
  OPTIONS$gaptol = gaptol
  OPTIONS$maxit = maxit
  OPTIONS$logsummary = 0
  OPTIONS$outputstats = 0
  OPTIONS$print = 0
  if (verbose)
    cat("Solving SDP ... ")
  sol = Rdsdp::dsdp(A, b, C, K, OPTIONS)
  if (verbose)
    cat("done. \n")
  if (!identical(sol$STATS$stype, "PDFeasible")) {
    warning("The SDP solver returned a non-feasible solution. Knockoffs may lose power.")
  }
  s = sol$y
  s[s < 0] = 0
  s[s > 1] = 1
  if (verbose)
    cat("Verifying that the solution is correct ... ")
  psd = 0
  s_eps = 1e-08
  while ((psd == 0) & (s_eps <= 0.1)) {
    if (is_posdef(2 * G - diag(s * (1 - s_eps), length(s)),
                  tol = 1e-09)) {
      psd = 1
    }
    else {
      s_eps = s_eps * 10
    }
  }
  s = s * (1 - s_eps)
  s[s < 0] = 0
  if (verbose)
    cat("done. \n")
  if (all(s == 0)) {
    warning("In creation of SDP knockoffs, procedure failed. Knockoffs will have no power.",
            immediate. = T)
  }
  return(s * diag(Sigma))
}

#' @keywords internal
create.solve_asdp <- function (Sigma, max.size = 500, gaptol = 1e-06, maxit = 1000,
                               verbose = FALSE)
{
  stopifnot(isSymmetric(Sigma))
  if (ncol(Sigma) <= max.size)
    return(create.solve_sdp(Sigma, gaptol = gaptol, maxit = maxit,
                            verbose = verbose))
  if (verbose)
    cat(sprintf("Dividing the problem into subproblems of size <= %s ... ",
                max.size))
  cluster_sol = divide.sdp(Sigma, max.size = max.size)
  n.blocks = max(cluster_sol$clusters)
  if (verbose)
    cat("done. \n")
  if (verbose)
    cat(sprintf("Solving %s smaller SDPs ... \n", n.blocks))
  s_asdp_list = list()
  if (verbose)
    pb <- utils::txtProgressBar(min = 0, max = n.blocks,
                                style = 3)
  for (k in 1:n.blocks) {
    s_asdp_list[[k]] = create.solve_sdp(as.matrix(cluster_sol$subSigma[[k]]),
                                        gaptol = gaptol, maxit = maxit)
    if (verbose)
      utils::setTxtProgressBar(pb, k)
  }
  if (verbose)
    cat("\n")
  p = dim(Sigma)[1]
  idx_count = rep(1, n.blocks)
  s_asdp = rep(0, p)
  for (j in 1:p) {
    cluster_j = cluster_sol$clusters[j]
    s_asdp[j] = s_asdp_list[[cluster_j]][idx_count[cluster_j]]
    idx_count[cluster_j] = idx_count[cluster_j] + 1
  }
  if (verbose)
    cat(sprintf("Combinining the solutions of the %s smaller SDPs ... ",
                n.blocks))
  tol = 1e-12
  maxitr = 1e+05
  gamma_range = seq(0, 1, len = 1000)
  options(warn = -1)
  gamma_opt = gtools::binsearch(function(i) {
    G = 2 * Sigma - gamma_range[i] * diag(s_asdp)
    lambda_min = RSpectra::eigs(G, 1, which = "SR", opts = list(retvec = FALSE,
                                                                maxitr = maxitr, tol = tol))$values
    if (length(lambda_min) == 0) {
      lambda_min = 1
    }
    lambda_min
  }, range = c(1, length(gamma_range)))
  s_asdp_scaled = gamma_range[min(gamma_opt$where)] * s_asdp
  options(warn = 0)
  if (verbose)
    cat("done. \n")
  if (verbose)
    cat("Verifying that the solution is correct ... ")
  if (!is_posdef(2 * Sigma - diag(s_asdp_scaled, length(s_asdp_scaled)))) {
    warning("In creation of approximate SDP knockoffs, procedure failed. Knockoffs will have no power.",
            immediate. = T)
    s_asdp_scaled = 0 * s_asdp_scaled
  }
  if (verbose)
    cat("done. \n")
  s_asdp_scaled
}

#'knockoff.filter
#'
#' @export
knockoff.filter <- function(X, y,
                            knockoffs=create.second_order,
                            statistic=stat.glmnet_coefdiff,
                            fdr=0.10,
                            offset=1
) {

  # Validate input types.
  if (is.data.frame(X)) {
    X.names = names(X)
    X = as.matrix(X, rownames.force = F)
  } else if (is.matrix(X)) {
    X.names = colnames(X)
  } else {
    stop('Input X must be a numeric matrix or data frame')
  }
  if (!is.numeric(X)) stop('Input X must be a numeric matrix or data frame')

  if (!is.factor(y) && !is.numeric(y)) {
    stop('Input y must be either of numeric or factor type')
  }
  if( is.numeric(y) ) y = as.vector(y)

  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }

  if (!is.function(knockoffs)) stop('Input knockoffs must be a function')
  if (!is.function(statistic)) stop('Input statistic must be a function')

  # Validate input dimensions
  n = nrow(X); p = ncol(X)
  stopifnot(length(y) == n)

  # If fixed-design knockoffs are being used, provive them with the response vector
  # in order to augment the data with new rows if necessary
  if( identical(knockoffs, create.fixed) )
    knockoffs = function(x) create.fixed(x, y=y)

  # Create knockoff variables
  knock_variables = knockoffs(X)

  # If fixed-design knockoffs are being used, update X and Y with the augmented observations (if present)
  if (is(knock_variables,"knockoff.variables")){
    X  = knock_variables$X
    Xk = knock_variables$Xk
    if(!is.null(knock_variables$y)) y  = knock_variables$y
    rm(knock_variables)
  } else if (is(knock_variables,"matrix")){
    Xk = knock_variables
    rm(knock_variables)
  } else {
    stop('Knockoff variables of incorrect type')
  }

  # Compute statistics
  W = statistic(X, Xk, y)

  # Run the knockoff filter
  t = knockoff.threshold(W, fdr=fdr, offset=offset)
  selected = sort(which(W >= t))
  if (!is.null(X.names))
    names(selected) = X.names[selected]

  # Package up the results.
  structure(list(call = match.call(),
                 X = X,
                 Xk = Xk,
                 y = y,
                 statistic = W,
                 threshold = t,
                 selected = selected),
            class = 'knockoff.result')
}


#' @keywords internal
knockoff.threshold <- function(W, fdr=0.10, offset=1) {
  if(offset!=1 && offset!=0) {
    stop('Input offset must be either 0 or 1')
  }
  ts = sort(c(0, abs(W)))
  ratio = sapply(ts, function(t)
    (offset + sum(W <= -t)) / max(1, sum(W >= t)))
  ok = which(ratio <= fdr)
  ifelse(length(ok) > 0, ts[ok[1]], Inf)
}

#' @keywords internal
create.second_order <- function (X, method = c("asdp", "equi", "sdp"), shrink = F)
{
  method = match.arg(method)
  mu = colMeans(X)
  if (!shrink) {
    Sigma = cov(X)
    if (!is_posdef(Sigma)) {
      shrink = TRUE
    }
  }
  if (shrink) {
    if (!requireNamespace("corpcor", quietly = T))
      stop("corpcor is not installed", call. = F)
    Sigma = tryCatch({
      suppressWarnings(matrix(as.numeric(corpcor::cov.shrink(X,
                                                             verbose = F)), nrow = ncol(X)))
    }, warning = function(w) {
    }, error = function(e) {
      stop("SVD failed in the shrinkage estimation of the covariance matrix. Try upgrading R to version >= 3.3.0")
    }, finally = {
    })
  }
  create.gaussian(X, mu, Sigma, method = method)
}
