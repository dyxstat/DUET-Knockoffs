#' FDR-Controlled Joint Feature Selection from Dual Microbiome data.
#'
#' @param W Count data.
#' @param class_K Vector data, representing different sources of count data.
#' @param data_x  Covariate matrix.
#' @param M The sequencing depth of count data.
#' @param y Response variable, typically binary.
#' @param T_var True correlation variables, default to NULL.
#' @param fdr The target FDR level, default 0.2.
#' @param offset value between 0 and 1.
#' @param test_statistic Type of single-platform test statistic, default 'DE', c('DE', 'GLM', 'RF')
#' @param filter_statistics Cross-platform aggregated knockoff contrast used to rank features:
#'   1 = product rule; 2 = max rule; 3 = sum rule.
#' @param test1 Test method, select when test_statistic = 'DE', default is 'wilcox.test', c('wilcox.test', 'ks.test').
#'
#' @return A list with:
#'   \item{test_stat}{Matrix of source‐specific test statistics.}
#'   \item{filter_stat}{Aggregated cross-platform statistic.}
#'   \item{S}{Selected feature indices.}
#'   \item{res}{Named vector of target FDR, empirical FDR, and power.}
#' @export
#'
#' @examples
#' library(DUETknockoff)
#' data(Zhu_count1)
#' data(Ye_count1)
#'
#' set.seed(42)
#'
#' n_group <- 2
#' n_data <- 2
#'
#' W1 <- Zhu_count1
#' W2 <- Ye_count1
#' W <- rbind(W1, W2)
#' M <- rowSums(W)
#' data_x <- NULL
#'
#' n_sam <- nrow(W1)
#' class_K <- factor(rep(1:n_data, each = n_sam))
#' y <- rep(c(rep(1, n_sam/2), rep(2, n_sam/2)), times = n_data)
#'
#' res.DUET.knockoff <- DUET_knockoff(W = W, M = M, class_K = class_K, y=y, data_x=NULL,
#'                                    T_var = NULL, test_statistic = "DE",
#'                                    filter_statistics=3, offset=1)
#' res.DUET.knockoff$test_stat[1:10]
#' res.DUET.knockoff$filter_stat[1:10]
DUET_knockoff <- function(W = W, class_K = NULL, data_x = NULL, M = NULL, y = y, T_var = NULL, fdr = 0.2, offset = 1,
                          test_statistic = "DE", filter_statistics = 3, test1 = "wilcox.test") {
  if (is.null(class_K)) {
    class_K <- rep(1, dim(W)[1])
  }
  if (is.null(data_x)) { # sum(data_x) == 1
    data_x <- as.data.frame(W[, c(1:3)])
    data_x[data_x != 0] <- 0
  }
  if (is.null(M)) {
    M <- apply(W, 1, sum)
  }
  name_data <- names(table(class_K))

  test_stat <- c()
  for (k1 in 1:length(name_data)) {
    sub <- class_K == name_data[k1]
    W_k <- W[sub, ]
    data_x_k <- data_x[sub, ]
    M_k <- M[sub]
    y_k <- y[sub]

    W_k <- apply(W_k, 2, function(col) {
      col_replace <- mean(col, na.rm = T)
      col[is.na(col)] <- col_replace
      if (any(is.infinite(col))) {
        col[is.infinite(col)] <- max(col[!is.infinite(col)])
      }
      if (any(is.na(col))) {
        col[is.na(col)] <- 0
      }
      return(col)
    })

    ## method = "ZINB"
    W_k_1 <- scDesign2_simulation(W_k, y_k)

    W_k_1 <- apply(W_k_1, 2, function(col) {
      col_replace <- mean(col, na.rm = T)
      col[is.na(col)] <- col_replace
      if (any(is.na(col))) {
        col[is.na(col)] <- 0
      }
      return(col)
    })

    if (test_statistic == "DE") {
      test_stat_k <- contrast_score_computation(W_k, W_k_1, y_k, test1)
    } else if (test_statistic == "GLM") {
      # test_statistic1 <- stat.glmnet_coefdiff
      random_m <- matrix(runif(dim(W_k)[1] * dim(W_k)[2], min = 0, max = 1), nrow = dim(W_k)[1])
      W_k_1 <- W_k_1 + random_m
      W_k <- W_k + random_m
      test_stat_k <- stat.glmnet_coefdiff(W_k, W_k_1, y_k)
    } else if (test_statistic == "RF") {
      # test_statistic1 <- stat.random_forest
      test_stat_k <- stat.random_forest(W_k, W_k_1, y_k)
    }
    test_stat <- rbind(test_stat, test_stat_k)
  }
  # print(test_stat)
  if (k1 == 1) {
    filter_stat <- test_stat
  } else {
    filter_stat <- switch(filter_statistics,
                    apply(test_stat, 2, cumprod)[dim(test_stat)[1], ],
                    apply(test_stat, 2, max),
                    apply(test_stat, 2, sum)
    )
  }

  result_fdr <- clipper_BC(filter_stat, fdr)
  S <- result_fdr$discovery

  tmp <- FDR_Power(S, T_var)
  res <- c(
    FDR_target = result_fdr$FDR,
    FDR_emp = tmp$FP[1],
    Power   = tmp$FP[2]
  )
  return(list(fdr_target = fdr, test_stat = test_stat, filter_stat = filter_stat, S = S, res = res))
}

