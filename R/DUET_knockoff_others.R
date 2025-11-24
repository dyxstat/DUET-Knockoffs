#' Importance statistics based on GLM with cross-validation.
#'
#' Computes importance statistics
#' \eqn{W_j = |Z_j| - |\tilde Z_j|} using cross‐validated
#' penalized generalized linear models from **glmnet**.
#'
#' @param X Numeric matrix (n × p). Original variables.
#' @param X_k Numeric matrix (n × p). Knockoff variables.
#' @param y Response vector of length n. Accepts numeric (gaussian, poisson),
#' factor/integer (binomial, multinomial), or survival matrix for \code{family="cox"}.
#' @param family response type, \code{"gaussian"}.
#' @param cores Number of cores used for parallel CV, default 2.

#' @param ... additional arguments specific to \code{glmnet}.
#' @return A vector of statistics \eqn{W} of length p.
#'
#' @details This function fits a penalized GLM using \code{glmnet} with 10-fold
#' cross-validation, extracts the fitted coefficients of
#' \code{(X, X_k)}, and constructs knockoff statistics by
#' comparing the magnitudes of corresponding coefficients.
#'
#' @family statistics
#'
#' @rdname stat.glmnet_coefdiff
#' @export
stat.glmnet_coefdiff <- function(X, X_k, y, family='gaussian', cores=2, ...) {

  if (!requireNamespace('glmnet', quietly=T))
    stop('glmnet is not installed', call.=F)
  parallel=T
  if (!requireNamespace('doParallel', quietly=T)) {
    warning('doParallel is not installed. Without parallelization, the statistics will be slower to compute', call.=F,immediate.=T)
    parallel=F
  }
  if (!requireNamespace('parallel', quietly=T)) {
    warning('parallel is not installed. Without parallelization, the statistics will be slower to compute.', call.=F,immediate.=T)
    parallel=F
  }

  # Register cores for parallel computation
  if (parallel) {
    ncores = parallel::detectCores(all.tests = TRUE, logical = TRUE)
    if( cores==2 ) {
      cores = min(2,ncores)
    }
    else {
      if (cores > ncores ) {
        warning(paste("The requested number of cores is not available. Using instead",ncores,"cores"),immediate.=T)
        cores = ncores
      }
    }
    if (cores>1) {
      doParallel::registerDoParallel(cores=cores)
      parallel = TRUE
    }
    else {
      parallel = FALSE
    }
  }

  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)

  p = ncol(X)

  # Compute statistics
  glmnet.coefs = cv_coeffs_glmnet(cbind(X.swap, Xk.swap), y, family=family, parallel=parallel, ...)
  if(family=="multinomial") {
    Z <- abs(glmnet.coefs[[1]][2:(2*p+1)])
    for(b in 2:length(glmnet.coefs)) {
      Z <- Z + abs(glmnet.coefs[[b]][2:(2*p+1)])
    }
  } else if (family=="cox") {
    Z <- glmnet.coefs[1:(2*p)]
  } else {
    Z <- glmnet.coefs[2:(2*p+1)]
  }
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig+p])

  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)

  # Stop the parallel cluster (if applicable)
  if (parallel) {
    if (cores>1) {
      doParallel::stopImplicitCluster()
    }
  }
  return(W)
}

#' @keywords internal
cv_coeffs_glmnet <- function(X, y, nlambda=500, intercept=T, parallel=T, ...) {
  # Standardize variables
  X = scale(X)

  n = nrow(X); p = ncol(X)

  if (!methods::hasArg(family) ) family = "gaussian"
  else family = list(...)$family

  if (!methods::hasArg(lambda) ) {
    if( identical(family, "gaussian") ) {
      if(!is.numeric(y)) {
        stop('Input y must be numeric.')
      }
      # Unless a lambda sequence is provided by the user, generate it
      lambda_max = max(abs(t(X) %*% y)) / n
      lambda_min = lambda_max / 2e3
      k = (0:(nlambda-1)) / nlambda
      lambda = lambda_max * (lambda_min/lambda_max)^k
    }
    else {
      lambda = NULL
    }
  }

  cv.glmnet.fit <- glmnet::cv.glmnet(X, y, lambda=lambda, intercept=intercept,
                                     standardize=F,standardize.response=F, parallel=parallel, ...)

  coef(cv.glmnet.fit, s = "lambda.min")
}

#' Importance statistics based on random forests
#'
#' Computes the difference statistic
#'   \deqn{W_j = |Z_j| - |\tilde{Z}_j|}
#' where \eqn{Z_j} and \eqn{\tilde{Z}_j} are the random forest feature importances
#' of the jth variable and its knockoff, respectively.
#'
#' @param X Numeric matrix (n × p). Original variables.
#' @param X_k Numeric matrix (n × p). Knockoff variables.
#' @param y vector of length n, containing the response variables. If a factor, classification is assumed,
#' otherwise regression is assumed.
#' @param ... additional arguments specific to \code{ranger} (see Details).
#'
#' @return A vector of statistics \eqn{W} of length p.
#' @details
#' This function calls \code{ranger} to compute impurity-based variable
#' importance for each column of \code{cbind(X, X_k)}. The knockoff statistic
#' compares the importance of each original variable with that of its knockoff
#' counterpart.
#'
#' For a complete list of the available additional arguments, see \code{\link[ranger]{ranger}}.
#'
#' @rdname stat.random_forest
#' @export
stat.random_forest <- function(X, X_k, y, ...) {
  if (!requireNamespace('ranger', quietly=T))
    stop('ranger is not installed', call.=F)

  # Randomly swap columns of X and Xk
  swap = rbinom(ncol(X),1,0.5)
  swap.M = matrix(swap,nrow=nrow(X),ncol=length(swap),byrow=TRUE)
  X.swap  = X * (1-swap.M) + X_k * swap.M
  Xk.swap = X * swap.M + X_k * (1-swap.M)

  # Compute statistics
  Z = random_forest_importance(cbind(X.swap, Xk.swap), y)
  p = ncol(X)
  orig = 1:p
  W = abs(Z[orig]) - abs(Z[orig+p])

  # Correct for swapping of columns of X and Xk
  W = W * (1-2*swap)
}

#' @keywords internal
random_forest_importance <- function(X, y, ...) {
  df = data.frame(y=y, X=X)
  rfFit = ranger::ranger(y~., data=df, importance="impurity", write.forest=F, ...)
  as.vector(rfFit$variable.importance)
}

#' glm.nb function: Estimation of the parameter (theta) for a Negative Binomial generalized linear model.
#' @keywords internal
glm.nb <- function (formula, data, weights, subset, na.action, start = NULL,
          etastart, mustart, control = glm.control(...), method = "glm.fit",
          model = TRUE, x = FALSE, y = TRUE, contrasts = NULL, ...,
          init.theta, link = log)
{
  loglik <- function(n, th, mu, y, w) sum(w * (lgamma(th +
                                                        y) - lgamma(th) - lgamma(y + 1) + th * log(th) + y *
                                                 log(mu + (y == 0)) - (th + y) * log(th + mu)))
  link <- substitute(link)
  fam0 <- if (missing(init.theta))
    do.call("poisson", list(link = link))
  else do.call("negative.binomial", list(theta = init.theta,
                                         link = link))
  mf <- Call <- match.call()
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  Terms <- attr(mf, "terms")
  if (method == "model.frame")
    return(mf)
  Y <- model.response(mf, "numeric")
  X <- if (!is.empty.model(Terms))
    model.matrix(Terms, mf, contrasts)
  else matrix(, NROW(Y), 0)
  w <- model.weights(mf)
  if (!length(w))
    w <- rep(1, nrow(mf))
  else if (any(w < 0))
    stop("negative weights not allowed")
  offset <- model.offset(mf)
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  n <- length(Y)
  if (!missing(method)) {
    if (!exists(method, mode = "function"))
      stop(gettextf("unimplemented method: %s", sQuote(method)),
           domain = NA)
    glm.fitter <- get(method)
  }
  else {
    method <- "glm.fit"
    glm.fitter <- stats::glm.fit
  }
  if (control$trace > 1)
    message("Initial fit:")
  fit <- glm.fitter(x = X, y = Y, weights = w, start = start,
                    etastart = etastart, mustart = mustart, offset = offset,
                    family = fam0, control = list(maxit = control$maxit,
                                                  epsilon = control$epsilon, trace = control$trace >
                                                    1), intercept = attr(Terms, "intercept") > 0)
  class(fit) <- c("glm", "lm")
  mu <- fit$fitted.values
  th <- as.vector(theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                           trace = control$trace > 2))
  if (control$trace > 1)
    message(gettextf("Initial value for 'theta': %f", signif(th)),
            domain = NA)
  fam <- do.call("negative.binomial", list(theta = th, link = link))
  iter <- 0
  d1 <- sqrt(2 * max(1, fit$df.residual))
  d2 <- del <- 1
  g <- fam$linkfun
  Lm <- loglik(n, th, mu, Y, w)
  Lm0 <- Lm + 2 * d1
  while ((iter <- iter + 1) <= control$maxit && (abs(Lm0 -
                                                     Lm)/d1 + abs(del)/d2) > control$epsilon) {
    eta <- g(mu)
    fit <- glm.fitter(x = X, y = Y, weights = w, etastart = eta,
                      offset = offset, family = fam, control = list(maxit = control$maxit,
                                                                    epsilon = control$epsilon, trace = control$trace >
                                                                      1), intercept = attr(Terms, "intercept") >
                        0)
    t0 <- th
    th <- theta.ml(Y, mu, sum(w), w, limit = control$maxit,
                   trace = control$trace > 2)
    fam <- do.call("negative.binomial", list(theta = th,
                                             link = link))
    mu <- fit$fitted.values
    del <- t0 - th
    Lm0 <- Lm
    Lm <- loglik(n, th, mu, Y, w)
    if (control$trace) {
      Ls <- loglik(n, th, Y, Y, w)
      Dev <- 2 * (Ls - Lm)
      message(sprintf("Theta(%d) = %f, 2(Ls - Lm) = %f",
                      iter, signif(th), signif(Dev)), domain = NA)
    }
  }
  if (!is.null(attr(th, "warn")))
    fit$th.warn <- attr(th, "warn")
  if (iter > control$maxit) {
    warning("alternation limit reached")
    fit$th.warn <- gettext("alternation limit reached")
  }
  if (length(offset) && attr(Terms, "intercept")) {
    null.deviance <- if (length(Terms))
      glm.fitter(X[, "(Intercept)", drop = FALSE], Y,
                 w, offset = offset, family = fam, control = list(maxit = control$maxit,
                                                                  epsilon = control$epsilon, trace = control$trace >
                                                                    1), intercept = TRUE)$deviance
    else fit$deviance
    fit$null.deviance <- null.deviance
  }
  class(fit) <- c("negbin", "glm", "lm")
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  Call$init.theta <- signif(as.vector(th), 10)
  Call$link <- link
  fit$call <- Call
  if (model)
    fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  if (x)
    fit$x <- X
  if (!y)
    fit$y <- NULL
  fit$theta <- as.vector(th)
  fit$SE.theta <- attr(th, "SE")
  fit$twologlik <- as.vector(2 * Lm)
  fit$aic <- -fit$twologlik + 2 * fit$rank + 2
  fit$contrasts <- attr(X, "contrasts")
  fit$xlevels <- .getXlevels(Terms, mf)
  fit$method <- method
  fit$control <- control
  fit$offset <- offset
  fit
}

#' theta.ml() function: Estimate (theta) of the Negative Binomial Distribution.
#' @keywords internal
theta.ml <- function (y, mu, n = sum(weights), weights, limit = 10, eps = .Machine$double.eps^0.25,
          trace = FALSE)
{
  score <- function(n, th, mu, y, w) sum(w * (digamma(th +
                                                        y) - digamma(th) + log(th) + 1 - log(th + mu) - (y +
                                                                                                           th)/(mu + th)))
  info <- function(n, th, mu, y, w) sum(w * (-trigamma(th +
                                                         y) + trigamma(th) - 1/th + 2/(mu + th) - (y + th)/(mu +
                                                                                                              th)^2))
  if (inherits(y, "lm")) {
    mu <- y$fitted.values
    y <- if (is.null(y$y))
      mu + residuals(y)
    else y$y
  }
  if (missing(weights))
    weights <- rep(1, length(y))
  t0 <- n/sum(weights * (y/mu - 1)^2)
  it <- 0
  del <- 1
  if (trace)
    message(sprintf("theta.ml: iter %d 'theta = %f'", it,
                    signif(t0)), domain = NA)
  while ((it <- it + 1) < limit && abs(del) > eps) {
    t0 <- abs(t0)
    del <- score(n, t0, mu, y, weights)/(i <- info(n, t0,
                                                   mu, y, weights))
    t0 <- t0 + del
    if (trace)
      message("theta.ml: iter", it, " theta =", signif(t0))
  }
  if (t0 < 0) {
    t0 <- 0
    warning("estimate truncated at zero")
    attr(t0, "warn") <- gettext("estimate truncated at zero")
  }
  if (it == limit) {
    warning("iteration limit reached")
    attr(t0, "warn") <- gettext("iteration limit reached")
  }
  attr(t0, "SE") <- sqrt(1/i)
  t0
}

#' mvrnorm function: Produces samples from the specified multivariate normal distribution
#' @keywords internal
mvrnorm <- function (n = 1, mu, Sigma, tol = 1e-06, empirical = FALSE,
          EISPACK = FALSE)
{
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p)))
    stop("incompatible arguments")
  if (EISPACK)
    stop("'EISPACK' is no longer supported by R", domain = NA)
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L])))
    stop("'Sigma' is not positive definite")
  X <- matrix(rnorm(p * n), n)
  if (empirical) {
    X <- scale(X, TRUE, FALSE)
    X <- X %*% svd(X, nu = 0)$v
    X <- scale(X, FALSE, TRUE)
  }
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*%
    t(X)
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma)))
    nm <- dn[[1L]]
  dimnames(X) <- list(nm, NULL)
  if (n == 1)
    drop(X)
  else t(X)
}

#' negative.binomial function: Fit a Negative Binomial generalized linear model with (theta).
#' @keywords internal
negative.binomial <- function (theta = stop("'theta' must be specified"), link = "log")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt"))
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"",
                       linktemp))
  }
  .Theta <- theta
  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", theta, envir = env)
  variance <- function(mu) mu + mu^2/.Theta
  validmu <- function(mu) all(mu > 0)
  dev.resids <- function(y, mu, wt) 2 * wt * (y * log(pmax(1,
                                                           y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))
  aic <- function(y, n, mu, wt, dev) {
    term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) +
      lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) -
      lgamma(.Theta + y)
    2 * sum(term * wt)
  }
  initialize <- expression({
    if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
    n <- rep(1, nobs)
    mustart <- y + (y == 0)/6
  })
  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    rnegbin(nsim * length(ftd), ftd, .Theta)
  }
  environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
  famname <- paste("Negative Binomial(", format(round(theta,
                                                      4)), ")", sep = "")
  structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                 validmu = validmu, valideta = stats$valideta, simulate = simfun),
            class = "family")
}


#' Coupla-based Count Simulation
#'
#' @param copula_result Coupla model parameters.
#' @param n             Number of simulated samples.
#' @param marginal      Marginal distribution: \code{"nb"} or \code{"Gamma"}.
#'
#' @return Simulated count matrix \code{p × n}
#' \code{copula_result}
#' @export
simulate_count_copula <- function(copula_result, n = 100,
                                  marginal = c('nb', 'Gamma')){
  marginal <- match.arg(marginal)

  p1 <- length(copula_result$taxon_sel1)
  if(p1 > 0){
    result1 <- mvrnorm(n = n, mu = rep(0.0, p1), Sigma = copula_result$cov_mat)
    result1 <- matrix(result1, nrow = n)
    result2 <- apply(result1, 2, pnorm)
    result2 <- matrix(result2, nrow = n)
  }
  p2 <- length(copula_result$taxon_sel2)
  if(marginal == 'nb'){
    if(p1 > 0){
      result31 <- t(sapply(1:p1, function(iter){
        param <- copula_result$marginal_param1[iter, ]
        qnbinom(pmax(0.0, result2[, iter] - param[1]) / (1-param[1]),
                size = param[2], mu = param[3])
      }))
    }
    if(p2 > 0){
      result32 <- t(sapply(1:p2, function(iter){
        param <- copula_result$marginal_param2[iter, ]
        rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
      }))
    }
  }else if(marginal == 'Gamma'){
    if(p1 > 0){
      result31 <- t(sapply(1:p1, function(iter){
        param <- copula_result$marginal_param1[iter, ]
        qgamma(max(0.0, result2[, iter] - param[1]), shape = param[2], scale = param[3] / param[2])
      }))
    }
    if(p2 > 0){
      result32 <- t(sapply(1:p2, function(iter){
        param <- copula_result$marginal_param2[iter, ]
        rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
      }))
    }
  }

  result <- matrix(0, nrow = p1 + p2 + length(copula_result$taxon_sel3), ncol = n)
  if(p1 > 0){
    result[copula_result$taxon_sel1, ] <- result31
  }
  if(p2 > 0){
    result[copula_result$taxon_sel2, ] <- result32
  }
  result
}



#' Simulate a count matrix for a single sam type based on a (w/o copula model)
#'
#' @param model_params A list containing marginal parameters for simulation
#' @inheritParams simulate_count_copula
#'
#' @return A \code{p × n} simulated count matrix based on marginal NB or Gamma
#' models, where \code{p} is derived from \code{model_params}.
#' @export
simulate_count_ind <- function(model_params, n = 100,
                               marginal = c('nb', 'Gamma')){
  marginal <- match.arg(marginal)

  if(model_params$sim_method == 'copula' || 'taxon_sel3' %in% names(model_params)){
    p1 <- length(model_params$taxon_sel1)
    p2 <- length(model_params$taxon_sel2)
    if(marginal == 'nb'){
      if(p1 > 0){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }
      if(p2 > 0){
        result32 <- t(sapply(1:p2, function(iter){
          param <- model_params$marginal_param2[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }
    }else if(marginal == 'Gamma'){
      if(p1 > 0){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
      if(p2 > 0){
        result32 <- t(sapply(1:p2, function(iter){
          param <- model_params$marginal_param2[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
    }
    result <- matrix(0, nrow = p1 + p2 + length(model_params$taxon_sel3), ncol = n)
    if(p1 > 0){
      result[model_params$taxon_sel1, ] <- result31
    }
    if(p2 > 0){
      result[model_params$taxon_sel2, ] <- result32
    }
  }else{
    p1 <- length(model_params$taxon_sel1)
    p2 <- length(model_params$taxon_sel2)
    result <- matrix(0, nrow = p1 + p2, ncol = n)
    if(p1 > 0){
      if(marginal == 'nb'){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rnbinom(n, size = param[2], mu = param[3])
        }))
      }else if(marginal == 'Gamma'){
        result31 <- t(sapply(1:p1, function(iter){
          param <- model_params$marginal_param1[iter, ]
          rbinom(n, 1, 1-param[1]) * rgamma(n, shape = param[2], scale = param[3] / param[2])
        }))
      }
      result[model_params$taxon_sel1, ] <- result31
    }
  }
  result
}

#' Simulate count matrix for experimental design
#'
#' @param model_params A list of fitted microbiome models (one per group or condition),
#'                     each containing parameters from either a copula-based or independent
#'                     marginal model.
#' @param total_count_new Total number of reads in the simulated count matrix.
#' @param n_sam_new      The total number of samples in the simulated count matrix.
#' @param sam_type_prop  The sam type proportion in the simulated count matrix.
#' @param total_count_old The total number of reads or in the original count matrix.
#' @param n_sam_old      The The total number of samples in the original count matrix.
#' @param sim_method      Specification of the type of model for data simulation.
#'                        Default value is 'copula', which selects the copula model.
#'                        'ind' will select the (w/o copula) model.
#' @param reseq_method    Strategy for adjusting sequencing depth:
#'                        \code{"mean_scale"} (default) rescales marginal means;
#'                        \code{"multinomial"} performs multinomial resampling to match
#'                        \code{total_count_new}.
#' @param sam_type     Logical, whether sams for each group should be sampled from a
#'                        multinomial distribution or follows the exact same proportion as
#'                        specified in \code{sam_type_prop}.
#' @return A matrix of shape p by n that contains the simulated count values. p is derived from
#' \code{model_params}.
#'
#' @export
simulate_count_scDesign2 <- function(model_params, n_sam_new, sam_type_prop = 1,
                                     total_count_new = NULL, total_count_old = NULL,
                                     n_sam_old = NULL, sim_method = c('copula', 'ind'),
                                     reseq_method = c('mean_scale', 'multinomial'),
                                     sam_type = FALSE){
  sim_method <- match.arg(sim_method)
  reseq_method <- match.arg(reseq_method)

  n_sam_vec <- sapply(model_params, function(x) x$n_sam)
  n_read_vec <- sapply(model_params, function(x) x$n_read)

  # if(is.null(total_count_new)) total_count_new <- sum(n_read_vec)
  # if(is.null(n_sam_new))      n_sam_new      <- sum(n_sam_vec)
  # if(is.null(sam_type_prop))  sam_type_prop  <- n_sam_vec
  if(is.null(total_count_old)) total_count_old <- sum(n_read_vec)
  if(is.null(n_sam_old))      n_sam_old      <- sum(n_sam_vec)

  if(is.null(total_count_new)) reseq_method <- 'mean_scale'

  if(length(model_params)!=length(sam_type_prop)){
    stop('sam type proportion should have the same length as the number of models.')
  }

  n_sam_type <- length(sam_type_prop)
  if(sam_type == TRUE){
    n_sam_each <- as.numeric(rmultinom(1, size = n_sam_new, prob = sam_type_prop))
  }else{
    sam_type_prop <- sam_type_prop / sum(sam_type_prop)
    n_sam_each <- round(sam_type_prop * n_sam_new)
    if(sum(n_sam_each) != n_sam_new){
      idx <- sample(n_sam_type, size = 1)
      n_sam_each[idx] <- n_sam_each[idx] + n_sam_new - sum(n_sam_each)
    }
  }

  p <- length(model_params[[1]]$taxon_sel1) + length(model_params[[1]]$taxon_sel2) +
    length(model_params[[1]]$taxon_sel3)
  new_count <- matrix(0, nrow = p, ncol = n_sam_new)
  if(reseq_method == 'mean_scale'){
    n_sam_each
    if(is.null(total_count_new)){
      r <- rep(1, n_sam_type)
    }else if(length(total_count_new) == 1){
      r <- rep(total_count_new / sum((total_count_old / n_sam_old) * n_sam_each),
               n_sam_type)
    }else{
      r <- (total_count_new / n_sam_new) / (total_count_old / n_sam_old)
    }
    for(iter in 1:n_sam_type)
      if(n_sam_each[iter] > 0){
        ulim <- sum(n_sam_each[1:iter])
        llim <- ulim - n_sam_each[iter] + 1
        params_new <- model_params[[iter]]
        params_new$marginal_param1[, 3] <- params_new$marginal_param1[, 3] * r[iter]
        if(sim_method == 'copula'){
          params_new$marginal_param2[, 3] <- params_new$marginal_param2[, 3] * r[iter]
          new_count[, llim:ulim] <- simulate_count_copula(params_new, n = n_sam_each[iter],
                                                          marginal = 'nb')
        }else if(sim_method == 'ind'){
          new_count[, llim:ulim] <- simulate_count_ind(params_new, n = n_sam_each[iter],
                                                       marginal = 'nb')
        }
      }
    if(is.null(names(model_params))){
      colnames(new_count) <- unlist(lapply(1:n_sam_type, function(x){rep(x, n_sam_each[x])}))
    }else{
      colnames(new_count) <- unlist(lapply(1:n_sam_type, function(x){
        rep(names(model_params)[x], n_sam_each[x])}))
    }
    return(new_count)
  }else if(reseq_method == 'multinomial'){
    for(iter in 1:n_sam_type){
      ulim <- sum(n_sam_each[1:iter])
      llim <- ulim - n_sam_each[iter] + 1
      if(sim_method == 'copula'){
        new_count[, llim:ulim] <- simulate_count_copula(model_params[[iter]],
                                                        n = n_sam_each[iter], marginal = 'Gamma')
      }else if(sim_method == 'ind'){
        new_count[, llim:ulim] <- simulate_count_ind(model_params[[iter]],
                                                     n = n_sam_each[iter], marginal = 'Gamma')
      }
    }

    new_count[which(is.infinite(new_count))] <- 0
    new_count[which(is.na(new_count))] <- 0

    bam_file <- sample(x = p*n_sam_new, size = total_count_new,
                       replace = TRUE, prob = as.vector(new_count))
    hist_result <- hist(bam_file, breaks = 0:(n_sam_new*p), plot = FALSE)
    result <- matrix(hist_result$counts, nrow = nrow(new_count))
    if(is.null(names(model_params))){
      colnames(result) <- unlist(lapply(1:n_sam_type, function(x){rep(x, n_sam_each[x])}))
    }else{
      colnames(result) <- unlist(lapply(1:n_sam_type, function(x){
        rep(names(model_params)[x], n_sam_each[x])}))
    }
    return(result)
  }
}


#' Fit the marginal distributions for each row of a count matrix
#'
#' @param x            Count matrix with \code{p × n}.
#' @param marginal     The types of marginal distribution.
#'                     Default value is 'auto_choose' which chooses between ZINB, NB, ZIP
#'                     and Poisson by a likelihood ratio test (lrt) and whether there is
#'                     underdispersion.
#'                     'zinb' will fit the ZINB model. If there is underdispersion, it
#'                     will choose between ZIP and Poisson by a lrt. Otherwise, it will try to
#'                     fit the ZINB model. If in this case, there is no zero at all or an error
#'                     occurs, it will fit an NB model instead.
#'                     'nb' fits the NB model that chooses between NB and Poisson depending
#'                     on whether there is underdispersion.
#'                     'poisson' simply fits the Poisson model.
#' @param pval_cutoff  Cutoff of p-value of the lrt that determines whether
#'                     there is zero inflation.
#' @param epsilon      Threshold value for preventing the transformed quantile
#'                     to collapse to 0 or 1.
#' @param jitter       Logical, whether a random projection should be performed in the
#'                     distributional transform.
#' @param DT           Logical, whether distributional transformed should be performed.
#'                     If set to FALSE, the returned object u will be NULL.
#'
#' @return             a list with the following components:
#'\describe{
#'  \item{params}{a matrix of shape p by 3. The values of each column are: the ZI proportion,
#'  the dispersion parameter (for Poisson, it's Inf), and the mean parameter.}
#'  \item{u}{NULL or a matrix of the same shape as x, which records the transformed quantiles,
#'  by DT.}
#'}
#' @export
fit_marginals <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          pval_cutoff = 0.05, epsilon = 1e-5,
                          jitter = TRUE, DT = TRUE){
  p <- nrow(x)
  n <- ncol(x)

  marginal <- match.arg(marginal)
  if(marginal == 'auto_choose'){
    params <- t(apply(x, 1, function(taxon){
      m <- mean(taxon)
      v <- var(taxon)
      if(m >= v){
        mle_Poisson <- glm(taxon ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- pscl::zeroinfl(taxon ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }else{
        mle_NB <- glm.nb(taxon ~ 1)
        if(min(taxon) > 0)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        else
          tryCatch({
            mle_ZINB <- pscl::zeroinfl(taxon ~ 1|1, dist = 'negbin')
            chisq_val <- 2 * (logLik(mle_ZINB) - logLik(mle_NB))
            pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
            if(pvalue < pval_cutoff)
              c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
            else
              c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          },
          error = function(cond){
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'zinb'){
    params <- t(apply(x, 1, function(taxon){
      m <- mean(taxon)
      v <- var(taxon)
      if(m >= v)
      {
        mle_Poisson <- glm(taxon ~ 1, family = poisson)
        tryCatch({
          mle_ZIP <- pscl::zeroinfl(taxon ~ 1|1, dist = 'poisson')
          chisq_val <- 2 * (logLik(mle_ZIP) - logLik(mle_Poisson))
          pvalue <- as.numeric(1 - pchisq(chisq_val, 1))
          if(pvalue < pval_cutoff)
            c(plogis(mle_ZIP$coefficients$zero), Inf, exp(mle_ZIP$coefficients$count))
          else
            c(0.0, Inf, m)
        },
        error = function(cond){
          c(0.0, Inf, m)})
      }
      else
      {
        if(min(taxon) > 0)
        {
          mle_NB <- glm.nb(taxon ~ 1)
          c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
        }
        else
          tryCatch({
            mle_ZINB <- pscl::zeroinfl(taxon ~ 1|1, dist = 'negbin')
            c(plogis(mle_ZINB$coefficients$zero), mle_ZINB$theta, exp(mle_ZINB$coefficients$count))
          },
          error = function(cond){
            mle_NB <- glm.nb(taxon ~ 1)
            c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
          })
      }
    }))
  }else if(marginal == 'nb'){
    params <- t(apply(x, 1, function(taxon){
      m <- mean(taxon)
      v <- var(taxon)
      if(m >= v){
        c(0.0, Inf, m)
      }else{
        mle_NB <- glm.nb(taxon ~ 1)
        c(0.0, mle_NB$theta, exp(mle_NB$coefficients))
      }
    }))
  }else if(marginal == 'poisson'){
    params <- t(apply(x, 1, function(taxon){
      c(0.0, Inf, mean(taxon))
    }))
  }

  if(DT){
    u <- t(sapply(1:p, function(iter){
      param <- params[iter, ]
      taxon <- unlist(x[iter,])
      prob0 <- param[1]
      u1 <- prob0 + (1 - prob0) * pnbinom(taxon, size = param[2], mu = param[3])
      u2 <- (prob0 + (1 - prob0) * pnbinom(taxon - 1, size = param[2], mu = param[3])) *
        as.integer(taxon > 0)
      if(jitter)
        v <- runif(n)
      else
        v <- rep(0.5, n)
      r <- u1 * v + u2 * (1 - v)
      idx_adjust <- which(1-r < epsilon)
      r[idx_adjust] <- r[idx_adjust] - epsilon
      idx_adjust <- which(r < epsilon)
      r[idx_adjust] <- r[idx_adjust] + epsilon

      r
    }))
  }else{
    u <- NULL
  }

  return(list(params = params, u = u))
}



#' Fit a Gaussian copula model for taxa counts
#'
#' @inheritParams fit_marginals
#' @param zp_cutoff Proportion threshold of zeros; taxa with zero proportion
#'   below this cutoff are included in the joint copula model.
#' @param min_nonzero_num Minimum number of non-zero counts required for a taxon
#'   to fit a marginal model.
#'
#' @return A list describing three taxon groups and the fitted copula model:
#' \describe{
#'   \item{cov_mat}{Correlation matrix of the Gaussian copula (group 1 taxa).}
#'   \item{marginal_param1}{Marginal parameters for group 1 taxa.}
#'   \item{marginal_param2}{Marginal parameters for group 2 taxa.}
#'   \item{taxa_sel1}{Indices of taxa included in the copula model.}
#'   \item{taxa_sel2}{Indices of taxa with only marginal models.}
#'   \item{taxa_sel3}{Indices of remaining taxa.}
#'   \item{zp_cutoff}{Input zero-proportion cutoff.}
#'   \item{min_nonzero_num}{Input minimum non-zero requirement.}
#'   \item{sim_method}{Character string \code{"copula"}.}
#' }
#' @export
fit_Gaussian_copula <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                                jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2){
  marginal <- match.arg(marginal)
  n <- ncol(x)
  p <- nrow(x)

  taxon_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })

  taxon_sel1 <- which(taxon_zero_prop < zp_cutoff)
  taxon_sel2 <- which(taxon_zero_prop < 1.0 - min_nonzero_num/n &
                       taxon_zero_prop >= zp_cutoff)
  taxon_sel3 <- (1:p)[-c(taxon_sel1, taxon_sel2)]

  if(length(taxon_sel1) > 0){
    marginal_result1 <- fit_marginals(x[taxon_sel1, , drop = FALSE], marginal, jitter = jitter, DT = TRUE)
    quantile_normal <- qnorm(marginal_result1$u)
    cov_mat <- cor(t(quantile_normal))
  }else{
    cov_mat = NULL
    marginal_result1 = NULL
  }

  if(length(taxon_sel2) > 0){
    marginal_result2 <- fit_marginals(x[taxon_sel2, , drop = FALSE], marginal, DT = FALSE)
  }else{
    marginal_result2 = NULL
  }
  return(list(cov_mat = cov_mat, marginal_param1 = marginal_result1$params,
              marginal_param2 = marginal_result2$params,
              taxon_sel1 = taxon_sel1, taxon_sel2 = taxon_sel2, taxon_sel3 = taxon_sel3,
              zp_cutoff = zp_cutoff, min_nonzero_num = min_nonzero_num,
              sim_method = 'copula', n_sam = n, n_read = sum(x)))
}



#' Fit an independent marginal model for taxa counts
#'
#' This function only fits the marginal distribution for each taxon.
#'
#' @inheritParams fit_marginals
#' @inheritParams fit_Gaussian_copula
#'
#' @return A list containing the fitted independent marginal model:
#' \describe{
#'   \item{marginal_param1}{Marginal parameters for taxa with at least
#'     \code{min_nonzero_num} non-zero counts.}
#'   \item{taxa_sel1}{Indices of taxa with sufficient non-zero counts.}
#'   \item{taxa_sel2}{Indices of taxa excluded from modeling.}
#'   \item{min_nonzero_num}{Input threshold for minimum non-zero counts.}
#'   \item{sim_method}{The string \code{"ind"} indicating an independent model.}
#' }
#' @export
fit_wo_copula <- function(x, marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                          jitter = TRUE, min_nonzero_num = 2){
  marginal <- match.arg(marginal)
  n <- ncol(x)
  p <- nrow(x)

  taxon_zero_prop <- apply(x, 1, function(y){
    sum(y < 1e-5) / n
  })

  taxon_sel1 <- which(taxon_zero_prop < 1.0 - min_nonzero_num/n)
  taxon_sel2 <- (1:p)[-taxon_sel1]

  if(length(taxon_sel1) > 0){
    marginal_result1 <- fit_marginals(x[taxon_sel1, ], marginal, jitter = jitter, DT = FALSE)
  }else{
    marginal_result1 = NULL
  }

  return(list(marginal_param1 = marginal_result1$params,
              taxon_sel1 = taxon_sel1, taxon_sel2 = taxon_sel2,
              min_nonzero_num = min_nonzero_num, sim_method = 'ind',
              n_sam = n, n_read = sum(x)))
}


#' Fit models for microbiome count matrix
#'
#' @param data_mat      A p × n count matrix (p taxa, n samples).
#'                      Column names should indicate sample-type (e.g., condition or group)
#'                      and must match \code{sam_type_sel}.
#' @param sam_type_sel  Character vector of selected sample-types for model fitting.
#' @param sim_method    Type of model to fit: 'copula' (default) or 'ind' (independent).
#' @inheritParams fit_Gaussian_copula
#' @param ncores         Number of parallel cores.
#' @return A list (length = number of sample types), where each element contains the fitted
#'         model for that sample type.
#'
#' @export
fit_model_scDesign2 <- function(data_mat, sam_type_sel, sim_method = c('copula', 'ind'),
                                marginal = c('auto_choose', 'zinb', 'nb', 'poisson'),
                                jitter = TRUE, zp_cutoff = 0.8,
                                min_nonzero_num = 2, ncores = 1){
  # if (!require(parallel)) install.packages("parallel")
  # library(parallel)
  sim_method <- match.arg(sim_method)
  marginal <- match.arg(marginal)

  if(sum(abs(data_mat - round(data_mat))) > 1e-5){
    warning('The entries in the input matrix are not integers. Rounding is performed.')
    data_mat <- round(data_mat)
  }

  if(sim_method == 'copula'){
    param <- parallel::mclapply(1:length(sam_type_sel), function(iter){
      fit_Gaussian_copula(data_mat[, colnames(data_mat) == sam_type_sel[iter]], marginal,
                          jitter = jitter, zp_cutoff = zp_cutoff,
                          min_nonzero_num = min_nonzero_num)
    }, mc.cores = ncores)
  }else if(sim_method == 'ind'){
    param <- parallel::mclapply(1:length(sam_type_sel), function(iter){
      fit_wo_copula(data_mat[, colnames(data_mat) == sam_type_sel[iter]], marginal,
                    jitter = jitter,
                    min_nonzero_num = min_nonzero_num)
    }, mc.cores = ncores)
  }

  names(param) <- sam_type_sel
  param
}




#' knockoff generation based on scDesign2 methods.
#'
#' @param W Count matrix.
#' @param class0 Response variable, binary variable.
#'
#' @return sim_count_copula: knockoff generation.
#' @export
#'
scDesign2_simulation <- function(W, class0) {
  # if (!require(scDesign2)) devtools::install_github("JSB-UCLA/scDesign2")
  # library(scDesign2)
  count_mat <- t(W) # W
  colnames(count_mat) <- class0
  sam_type_sel <- unique(class0)
  # submat <- get_submat(count_mat, sam_type_sel)
  n_sam_new <- ncol(count_mat)
  sam_type_prop <- table(colnames(count_mat))[sam_type_sel]
  copula_result <- fit_model_scDesign2(count_mat, sam_type_sel,
                                       sim_method = "copula", marginal = "zinb",
                                       ncores = 1
  ) # length(sam_type_sel)
  sim_count_copula <- simulate_count_scDesign2(copula_result, n_sam_new,
                                               sim_method = "copula",
                                               # marginal = 'zinb',
                                              sam_type_prop = sam_type_prop
  )
  return(t(sim_count_copula))
}


#' Calculation of test statistics based on DE (wilcoxon, ks.test).
#'
#' @param W Count matrix;
#' @param result3 knockoff generation.
#' @param class0 Response variable, binary variable or others.
#' @param test1 Selection of test statistics, c("wilcox.test", "ks.test")
#'
#' @return contrast_score: contrast score.
#' @export
contrast_score_computation <- function(W, result3, class0, test1) {
  # if (!require(uwot)) install.packages("uwot")
  # if (!require(cluster)) install.packages("cluster")
  # library(uwot)
  # library(cluster)

  # result3 <- result_konckoffs$Xk

  test2 <- c("wilcox.test", "ks.test")
  test1 <- match(test1, test2)

  allmean <- mean(result3, na.rm = TRUE)
  result3_1 <- apply(result3, 2, function(col) {
    # col_mean <- mean(col, na.rm = TRUE)
    col_replace <- sample(c(0, allmean), sum(is.na(col)), replace = TRUE, prob = c(0.8, 0.2))
    col[is.na(col)] <- col_replace
    return(col)
  })

  umap_result <- uwot::umap(result3_1, n_components = 3)
  # kmeans_result <- kmeans(umap_result$layout, centers = 2)
  kmeans_result <- stats::kmeans(umap_result, centers = 2)
  class1 <- kmeans_result$cluster

  contrast_score <- sapply(1:dim(result3)[2], function(i_wi) {
    # odd_vector <- result3[class1==1,i_Wi]
    # even_vector <- result3[class1==2,i_Wi]
    # KS.test(odd_vector,even_vector)$p.value
    P <- sort(result3[class1 == 1, i_wi])
    Q <- sort(result3[class1 == 2, i_wi])

    P1 <- sort(W[class0 == 1, i_wi])
    Q1 <- sort(W[class0 == 2, i_wi])

    if (sum(P) + sum(Q) == 0 | sum(P1) + sum(Q1) == 0) {
      result_p <- 0
    } else {
      result_p <- switch(test1,
                         log10(stats::wilcox.test(P, Q)$p.value / stats::wilcox.test(P1, Q1)$p.value),
                         log10(stats::ks.test(P, Q)$p.value / stats::ks.test(P1, Q1)$p.value)
      )
    }
  })
  return(contrast_score)
}

#' FDR control method based on clipper.
#'
#' @param contrastScore Contrast Score or fliter statistics;
#' @param FDR Scalar or vector, target control level.
#'
#' @return
#'  FDR: If is vector, output FDR value of result;
#'  thre: Control thresholds for contrast scores;
#'  q: The converted values of the scores were compared;
#'  discovery: Results of variable selection at the FDR level.
#' @export
#'
clipper_BC <- function(contrastScore, FDR) {
  # contrastScore <- contrast_score
  contrastScore[is.na(contrastScore)] <- 0 # impute missing contrast scores with 0
  c_abs <- abs(contrastScore[contrastScore != 0])
  c_abs <- sort(unique(c_abs))

  i <- 1
  emp_fdp <- rep(NA, length(c_abs))
  emp_fdp[1] <- 1
  while (i <= length(c_abs)) {
    # print(i)
    t <- c_abs[i]
    emp_fdp[i] <- min((1 + sum(contrastScore <= -t)) / sum(contrastScore >= t), 1)
    if (i >= 2) {
      emp_fdp[i] <- min(emp_fdp[i], emp_fdp[i - 1])
    }
    i <- i + 1
  }

  c_abs <- c_abs[!is.na(emp_fdp)]
  emp_fdp <- emp_fdp[!is.na(emp_fdp)]
  q <- emp_fdp[match(contrastScore, c_abs)]
  q[which(is.na(q))] <- 1

  for (FDR_i in FDR) {
    re_i <- list(
      FDR = FDR_i,
      # thre = thre,
      q = q,
      discovery = NULL
    )
    if (sum(emp_fdp <= FDR_i) > 0) {
      thre <- c_abs[min(which(emp_fdp <= FDR_i))]
      re_i <- list(
        FDR = FDR_i,
        thre = thre,
        q = q,
        discovery = which(contrastScore >= thre)
      )
      # re_i = c(FDR_i,thre,round(which(contrastScore >= thre)))
      break
    }
  }
  return(re_i)

}

#' Calculation of FDR and Power
#'
#' @param a_1 Results of selected variables;
#' @param b_1 True significant variable.
#'
#' @return FDR and Power.
#' @export
#'
FDR_Power <- function(a_1, b_1 = NULL) {
  if (is.null(b_1)) {
    return(list(FP = c(-1, -1)))
  } else {
    if (length(b_1) == 1) {
      b_1 <- 1:b_1
    }
    fdr_moni_k <- sum(!(a_1 %in% b_1)) / max(length(a_1), 1)
    power_moni_k <- sum(b_1 %in% a_1) / length(b_1)
    return(list(FP = c(fdr_moni_k, power_moni_k)))
  }
}
