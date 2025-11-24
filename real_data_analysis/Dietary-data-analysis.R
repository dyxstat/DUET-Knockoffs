

library(readxl)
library(writexl)
library(ANCOMBC)
library(LOCOM)
library(permute)
library(ZIPGSK)

set.seed(1003)
seed <- 1003

fdr_target <- 0.1
n_rej_stop <- 100
n_cores <- 4

# W <- read.csv("C:/Users/AW/Desktop/real-data ZINBSK/W.csv", row.names = 1)
# M <- read.csv("C:/Users/AW/Desktop/real-data ZINBSK/M.csv", header = FALSE)
# y <- read.csv("C:/Users/AW/Desktop/real-data ZINBSK/y.csv", header = FALSE)
# data_x <- read.csv("C:/Users/AW/Desktop/real-data ZINBSK/data_x.csv")
# class_K <- read.csv("C:/Users/AW/Desktop/real-data ZINBSK/class_K.csv", header = FALSE)


W <- read.csv("C:/Users/AW/Desktop/real-data ZINBSK/Diet/W.csv", row.names = 1)
M <- read.csv("C:/Users/AW/Desktop/real-data ZINBSK/Diet/M.csv", header = FALSE)
y <- read.csv("C:/Users/AW/Desktop/real-data ZINBSK/Diet/y.csv", header = FALSE)
data_x <- read.csv("C:/Users/AW/Desktop/real-data ZINBSK/Diet/data_x.csv")
class_K <- read.csv("C:/Users/AW/Desktop/real-data ZINBSK/Diet/class_K.csv", header = FALSE)

# ==== Com2seq / Locom ==== #
W1 <- as.matrix(W[1:76, ])
W2 <- as.matrix(W[77:152, ])
rownames(W2) <- sub("1$", "", rownames(W2))
W_pool <- W1 + W2


Y1 <- as.numeric(y[1:76, 1])
Y2 <- as.numeric(y[77:152, 1])

Y1 <- ifelse(Y1 == 2, 1, 0)
Y2 <- ifelse(Y2 == 2, 1, 0)

C1 <- as.matrix(data_x[1:76, ])
C2 <- as.matrix(data_x[77:152, ])


# locom_16s, locom_sms, locom_pool
res.locom1 <- locom(W1, Y=Y1, C=C1, n.perm.max = 1000,
                    fdr.nominal = fdr_target, seed = seed,
                    n.cores = n_cores, n.rej.stop = n_rej_stop)

res.locom2 <- locom(W2, Y=Y2, C=C1, n.perm.max = 1000,
                    fdr.nominal = fdr_target, seed = seed,
                    n.cores = n_cores, n.rej.stop = n_rej_stop)

res.locom.pool <- locom(W_pool, Y=Y1, C=C1, n.perm.max = 1000,
                        fdr.nominal = fdr_target, seed = seed,
                        n.cores = n_cores, n.rej.stop = n_rej_stop)

res.Com2seq <- Com2seq(
  table1 = W1, table2 = W2, Y1=Y1, Y2=Y2, C1=C1, C2=C2,
  fdr.nominal = fdr_target, n.cores = n_cores,
  n.perm.max = 1000, n.rej.stop = n_rej_stop,
  filter.thresh = 0, seed = seed
)


# ==== ANCOM-BC ==== #
meta_W <- data.frame(
  Group=factor(Y1), Cov1=C1[,1], Cov2=C1[,2], Cov3=C1[,3], Cov4=C1[,4],
  row.names=rownames(W1)
)

res.ancombc1 <- ancombc2(
  data = t(W1), meta_data = meta_W, taxa_are_rows = TRUE,
  fix_formula = "Group + Cov1 + Cov2 + Cov3",
  p_adj_method = "holm", pseudo_sens = F,
  prv_cut = 0.20, group = "Group",
  alpha = fdr_target, n_cl = 1,
  global = TRUE, verbose = TRUE
)

res.ancombc2 <- ancombc2(
  data = t(W2), meta_data = meta_W, taxa_are_rows = TRUE,
  fix_formula = "Group + Cov1 + Cov2 + Cov3",
  p_adj_method = "holm", pseudo_sens = T,
  prv_cut = 0.20, group = "Group",
  alpha = fdr_target, n_cl = n_cores,
  global = TRUE, verbose = TRUE
)

res.ancombc.pool <- ancombc2(
  data = t(W_pool), meta_data = meta_W, taxa_are_rows = TRUE,
  fix_formula = "Group + Cov1 + Cov2 + Cov3",
  p_adj_method = "holm", pseudo_sens = T,
  prv_cut = 0.20, group = "Group",
  alpha = fdr_target, n_cl = n_cores,
  global = TRUE, verbose = TRUE
)

# ==== Wilcoxon ==== #
group0_idx <- which(Y1 == 0)
group1_idx <- which(Y1 == 1)

p_wilcox1 <- apply(W1, 2, function(x) {
  tryCatch({
    wilcox.test(
      x[group0_idx], x[group1_idx],
      alternative = "two.sided",
      exact = FALSE
    )$p.value
  }, error = function(e) {
    return(NA)  # 如果检验失败，返回 NA
  })
})

p_wilcox2 <- apply(W2, 2, function(x) {
  tryCatch({
    wilcox.test(
      x[group0_idx], x[group1_idx],
      alternative = "two.sided",
      exact = FALSE
    )$p.value
  }, error = function(e) {
    return(NA)
  })
})

p_wilcox.pool <- apply(W_pool, 2, function(x) {
  tryCatch({
    wilcox.test(
      x[group0_idx], x[group1_idx],
      alternative = "two.sided",
      exact = FALSE
    )$p.value
  }, error = function(e) {
    return(NA)
  })
})

# ==== Cauchy p-value combination ==== #
p.locom1 <- res.locom1$p.otu  # LOCOM
p.locom2 <- res.locom2$p.otu

p.ancombc1 <- res.ancombc1$res$p_Group1 # ANCOMBC
names(p.ancombc1) <- res.ancombc1$res$taxon

p.ancombc2 <- res.ancombc2$res$p_Group1  # ANCOMBC
names(p.ancombc2) <- res.ancombc2$res$taxon

p.wilcox1 <- p_wilcox1
p.wilcox2 <- p_wilcox2

# ==== Loop for Cauchy combination ==== #

methods <- c("locom", "ancombc", "wilcox")
results <- list()

for(method in methods){
  cat("\n=== Processing", toupper(method), "Cauchy Combination ===\n")
  # obtain p-values
  p1 <- get(paste0("p.", method, "1"))
  p2 <- get(paste0("p.", method, "2"))

  if (is.matrix(p1)) p1 <- as.vector(p1)
  if (is.matrix(p2)) p2 <- as.vector(p2)

  # Cauchy combination
  name1 <- names(p1)
  name2 <- names(p2)
  common.name <- intersect(name1, name2)

  if (length(common.name) == 0) {
    cat("  Warning: No common taxa between datasets for", method, "\n")

    # if no common taxa, combine all the names
    p.comp <- c(p1, p2)
    p.comp.name <- c(name1, name2)
  } else {
    j.mat1 <- match(common.name, name1)
    j.mat2 <- match(common.name, name2)

    p.both <- rbind(p1[j.mat1], p2[j.mat2])
    p.comp <- pcauchy(apply(tan((0.5 - p.both)*pi), 2, mean),
                      lower.tail = FALSE)
    p.comp.name <- common.name

    # add non-intersect taxa
    if (length(j.mat1) < length(p1)) {
      p.comp <- c(p.comp, p1[-j.mat1])
      p.comp.name <- c(p.comp.name, name1[-j.mat1])
    }
    if (length(j.mat2) < length(p2)) {
      p.comp <- c(p.comp, p2[-j.mat2])
      p.comp.name <- c(p.comp.name, name2[-j.mat2])
    }
  }

  # BH method
  p.comp <- matrix(p.comp, nrow=1)
  q.comp <- matrix(p.adjust(p.comp, method = "BH"), nrow=1)
  colnames(p.comp) <- p.comp.name
  colnames(q.comp) <- p.comp.name

  # Gloabl p-value
  if (ncol(p.comp) > 0) {
    p.global <- pcauchy(mean(tan((0.5 - p.comp)*pi)), lower.tail = FALSE)
  } else {
    p.global <- NA
  }

  # Extract DA taxa
  cauchy_da <- colnames(q.comp)[which(q.comp < fdr_target)]

  # Extract corresponding q-values for DA taxa
  cauchy_q <- if(length(cauchy_da) > 0) {
    as.numeric(q.comp[, cauchy_da])
  } else {
    numeric(0)
  }

  # output
  results[[method]] <- list(
    p.comp = p.comp,
    q.comp = q.comp,
    p.global = p.global,
    cauchy_da = cauchy_da,
    cauchy_q = cauchy_q
  )

  cat(method, "Cauchy DA taxa:", length(cauchy_da), "\n")
  cat(method, "Global p-value:", p.global, "\n")
}

# Extract output and q-value for Cauchy Combination
DA_Cauchy_Locom <- results$locom$cauchy_da
Q_Cauchy_Locom <- results$locom$cauchy_q

DA_Cauchy_ancombc <- results$ancombc$cauchy_da
Q_Cauchy_ancombc <- results$ancombc$cauchy_q

DA_Cauchy_wilcox <- results$wilcox$cauchy_da
Q_Cauchy_wilcox <- results$wilcox$cauchy_q



# ==== Extract pooled DA and Com2seq taxa==== #
# Com2seq_omni
DA_Com2seq <- res.Com2seq$detected.taxa.omni

Q_Com2seq <- if(length(DA_Com2seq) > 0 && !is.null(res.Com2seq$q.taxa.omni)) {
  as.numeric(res.Com2seq$q.taxa.omni[, DA_Com2seq])
} else {
  numeric(0)
}
# Locom_16s
DA_Locom_16s <- res.locom1$detected.otu
Q_Locom_16s <- if(length(DA_Locom_16s) > 0) {
  as.numeric(res.locom1$q.otu[, DA_Locom_16s])
} else {
  numeric(0)
}

# Locom_shotgun
DA_Locom_shotgun <- res.locom2$detected.otu
Q_Locom_shotgun <- if(length(DA_Locom_shotgun) > 0) {
  as.numeric(res.locom2$q.otu[, DA_Locom_shotgun])
} else {
  numeric(0)
}

# Locom_pooled
DA_Pool_Locom <- res.locom.pool$detected.otu
Q_Pool_Locom <- if(length(DA_Pool_Locom) > 0) {
  as.numeric(res.locom.pool$q.otu[, DA_Pool_Locom])
} else {
  numeric(0)
}

# ANCOMBC_pooled
DA_Pool_ancombc <- res.ancombc.pool$res$diff_abn$taxon[
  res.ancombc.pool$res$diff_abn$Group1 == TRUE]

Q_Pool_ancombc <- if(length(DA_Pool_ancombc) > 0) {
  q_all <- res.ancombc.pool$res$q_Group1
  names(q_all) <- res.ancombc.pool$res$taxon
  as.numeric(q_all[DA_Pool_ancombc])
} else {
  numeric(0)
}

# Wilcoxon_pooled
q_wilcox.pool <- p.adjust(p_wilcox.pool, method = "BH")
DA_Pool_wilcox <- names(q_wilcox.pool)[
  !is.na(q_wilcox.pool) & q_wilcox.pool < fdr_target
]

Q_wilcox_pool <- if(length(DA_Pool_wilcox) > 0) {
  as.numeric(q_wilcox.pool[DA_Pool_wilcox])
} else {
  numeric(0)
}

# ==== ZINBSK ==== #
# M <- as.numeric(M[,1])
# class_K <- as.numeric(class_K[,1])
# y <- as.numeric(y[,1])


res.ZINBSK <- ZIPG_SK_other(W = W, M = M, class_K = class_K, data_x = data_x,
                            fdr = 0.1, method = "ZINB", y = y, T_var = NULL,
                            test_statistic = "GLM", filter_statistics = 2)
zinbsk_idx <- res.ZINBSK$S
DA_zinbsk <- colnames(W)[zinbsk_idx]
DA_zinbsk

W_zinbsk <- if(length(zinbsk_idx) > 0){
  res.ZINBSK$c_w_b[zinbsk_idx]
} else {
  numeric(0)
}

# ==== Summary Output ==== #

Taxa_list <- list(
  ZINBSK = DA_zinbsk,
  Com2seq_omni = DA_Com2seq,
  LOCOM_16s = DA_Locom_16s,
  LOCOM_shotgun = DA_Locom_shotgun,
  LOCOM_Com_p = DA_Cauchy_Locom,
  LOCOM_Com_count = DA_Pool_Locom,
  ANCOMBC_Com_p = DA_Cauchy_ancombc,
  ANCOMBC_Com_count = DA_Pool_ancombc,
  Wilcoxon_Com_p = DA_Cauchy_wilcox,
  Wilcoxon_Com_count = DA_Pool_wilcox
)

Q_value_list <- list(
  ZINBSK = W_zinbsk,
  Com2seq_omni = Q_Com2seq,
  LOCOM_16s = Q_Locom_16s,
  LOCOM_shotgun = Q_Locom_shotgun,
  LOCOM_Com_p = Q_Cauchy_Locom,
  LOCOM_Com_count = Q_Pool_Locom,
  ANCOMBC_Com_p = Q_Cauchy_ancombc,
  ANCOMBC_Com_count = Q_Pool_ancombc,
  Wilcoxon_Com_p = Q_Cauchy_wilcox,
  Wilcoxon_Com_count = Q_wilcox_pool
)

# ==== Create Wide Format Table ==== #
# all_taxa <- unique(unlist(Taxa_list))
# methods_names <- names(Taxa_list)
# q_value_matrix <- matrix(NA, nrow = length(methods_names),
#                          ncol = length(all_taxa))
# rownames(q_value_matrix) <- methods_names
# colnames(q_value_matrix) <- all_taxa
#
# for (method in methods_names) {
#   taxa <- Taxa_list[[method]]
#   q_vals <- Q_value_list[[method]]
#
#   if (length(taxa) > 0 && length(q_vals) > 0) {
#     for (j in 1:length(taxa)) {
#       taxon <- taxa[j]
#       q_value_matrix[method, taxon] <- q_vals[j]
#     }
#   }
# }
#
# DA_table <- as.data.frame(q_value_matrix)
# DA_table <- cbind(Method = rownames(DA_table), DA_table)
#
# write_xlsx(DA_table, "DA_Q_values_Summary.xlsx")
# cat("\n=== DA Summary Table Preview ===\n")
# print(DA_table[, 1:min(6, ncol(DA_table))])

