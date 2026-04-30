#############################################################
# Dietary niche — Single-run detected taxa
# Methods: DUETKnockoffs, Com2seq,
#          Locom_Com_P, Locom_Com_Count, Locom_16s, Locom_shotgun
#############################################################

# Set working directory to script location (RStudio only)
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(LOCOM)
library(DUETKnockoffs)

seed <- 2024
fdr_target <- 0.1
n_rej_stop <- 100

# ==== Load data ==== #
W       <- read.csv("Qiita_ID_11212/W.csv", row.names = 1)
M       <- read.csv("Qiita_ID_11212/M.csv", header = FALSE)
y       <- read.csv("Qiita_ID_11212/y.csv", header = FALSE)
data_x  <- read.csv("Qiita_ID_11212/data_x.csv")
class_K <- read.csv("Qiita_ID_11212/class_K.csv", header = FALSE)

set.seed(seed)

# ==== Split 16S / Shotgun / Pool ==== #
W1 <- as.matrix(W[1:76, ])
W2 <- as.matrix(W[77:152, ])
rownames(W2) <- rownames(W1)
W_pool <- W1 + W2

Y1 <- as.numeric(y[1:76, 1])
Y2 <- as.numeric(y[77:152, 1])
Y1 <- ifelse(Y1 == 2, 1, 0)
Y2 <- ifelse(Y2 == 2, 1, 0)

C1 <- as.matrix(data_x[1:76, ])
C2 <- as.matrix(data_x[77:152, ])

# ==== LOCOM: 16S / shotgun / pool ==== #
res.locom1 <- locom(W1, Y = Y1, C = C1, n.perm.max = 1000,
                    fdr.nominal = fdr_target, seed = seed,
                    n.cores = 1, n.rej.stop = n_rej_stop)

res.locom2 <- locom(W2, Y = Y2, C = C1, n.perm.max = 1000,
                    fdr.nominal = fdr_target, seed = seed,
                    n.cores = 1, n.rej.stop = n_rej_stop)

res.locom.pool <- locom(W_pool, Y = Y1, C = C1, n.perm.max = 1000,
                        fdr.nominal = fdr_target, seed = seed,
                        n.cores = 1, n.rej.stop = n_rej_stop)

# ==== Com2seq ==== #
res.Com2seq <- Com2seq(
  table1 = W1, table2 = W2, Y1 = Y1, Y2 = Y2, C1 = C1, C2 = C2,
  fdr.nominal = fdr_target, n.cores = 1,
  n.perm.max = 1000, n.rej.stop = n_rej_stop,
  filter.thresh = 0, seed = seed
)

# ==== Cauchy combination for LOCOM (Locom_Com_P) ==== #
p.locom1 <- res.locom1$p.otu
p.locom2 <- res.locom2$p.otu

if (is.matrix(p.locom1)) p.locom1 <- as.vector(p.locom1)
if (is.matrix(p.locom2)) p.locom2 <- as.vector(p.locom2)

name1 <- names(p.locom1)
name2 <- names(p.locom2)
common.name <- intersect(name1, name2)

if (length(common.name) == 0) {
  cat("  Warning: No common taxa between 16S and shotgun for LOCOM\n")
  p.comp <- c(p.locom1, p.locom2)
  p.comp.name <- c(name1, name2)
} else {
  j.mat1 <- match(common.name, name1)
  j.mat2 <- match(common.name, name2)

  p.both <- rbind(p.locom1[j.mat1], p.locom2[j.mat2])
  p.comp <- pcauchy(apply(tan((0.5 - p.both) * pi), 2, mean),
                    lower.tail = FALSE)
  p.comp.name <- common.name

  if (length(j.mat1) < length(p.locom1)) {
    p.comp <- c(p.comp, p.locom1[-j.mat1])
    p.comp.name <- c(p.comp.name, name1[-j.mat1])
  }
  if (length(j.mat2) < length(p.locom2)) {
    p.comp <- c(p.comp, p.locom2[-j.mat2])
    p.comp.name <- c(p.comp.name, name2[-j.mat2])
  }
}

p.comp <- matrix(p.comp, nrow = 1)
q.comp <- matrix(p.adjust(p.comp, method = "BH"), nrow = 1)
colnames(p.comp) <- p.comp.name
colnames(q.comp) <- p.comp.name

if (ncol(p.comp) > 0) {
  p.global.locom <- pcauchy(mean(tan((0.5 - p.comp) * pi)), lower.tail = FALSE)
} else {
  p.global.locom <- NA
}

DA_Cauchy_Locom <- colnames(q.comp)[which(q.comp < fdr_target)]
Q_Cauchy_Locom  <- if (length(DA_Cauchy_Locom) > 0) {
  as.numeric(q.comp[, DA_Cauchy_Locom])
} else numeric(0)

cat("LOCOM Cauchy DA taxa:", length(DA_Cauchy_Locom), "\n")
cat("LOCOM Global p-value:", p.global.locom, "\n")

# ==== Extract DA taxa per method ==== #

# Com2seq_omni
DA_Com2seq <- res.Com2seq$detected.taxa.omni
Q_Com2seq  <- if (length(DA_Com2seq) > 0 && !is.null(res.Com2seq$q.taxa.omni)) {
  as.numeric(res.Com2seq$q.taxa.omni[, DA_Com2seq])
} else numeric(0)

# Locom_16s
DA_Locom_16s <- res.locom1$detected.otu
Q_Locom_16s  <- if (length(DA_Locom_16s) > 0) {
  as.numeric(res.locom1$q.otu[, DA_Locom_16s])
} else numeric(0)

# Locom_shotgun
DA_Locom_shotgun <- res.locom2$detected.otu
Q_Locom_shotgun  <- if (length(DA_Locom_shotgun) > 0) {
  as.numeric(res.locom2$q.otu[, DA_Locom_shotgun])
} else numeric(0)

# Locom_Com_count (pooled counts)
DA_Pool_Locom <- res.locom.pool$detected.otu
Q_Pool_Locom  <- if (length(DA_Pool_Locom) > 0) {
  as.numeric(res.locom.pool$q.otu[, DA_Pool_Locom])
} else numeric(0)

# ==== DUET-Knockoffs ==== #
res.DUETKnockoffs <- DUET_Knockoffs(
  W = W, M = M, class_K = class_K, data_x = data_x,
  fdr = fdr_target, y = y, T_var = NULL,
  test_statistic = "DE", filter_statistics = 2
)
knockoff_idx <- res.DUETKnockoffs$S
DA_knockoff  <- colnames(W)[knockoff_idx]
W_knockoff   <- if (length(knockoff_idx) > 0) {
  res.DUETKnockoffs$filter_stat[knockoff_idx]
} else numeric(0)

cat("DUET-Knockoffs DA taxa:", length(DA_knockoff), "\n")

# ==== Summary Output ==== #
Taxa_list <- list(
  `DUET-Knockoffs`    = DA_knockoff,
  Com2seq         = DA_Com2seq,
  Locom_Com_P     = DA_Cauchy_Locom,
  Locom_Com_Count = DA_Pool_Locom,
  Locom_16s       = DA_Locom_16s,
  Locom_shotgun   = DA_Locom_shotgun
)

Q_value_list <- list(
  `DUET-Knockoffs`    = W_knockoff,
  Com2seq         = Q_Com2seq,
  Locom_Com_P     = Q_Cauchy_Locom,
  Locom_Com_Count = Q_Pool_Locom,
  Locom_16s       = Q_Locom_16s,
  Locom_shotgun   = Q_Locom_shotgun
)

cat(sprintf(
  "seed = %d | DUET=%d | Com2seq=%d | LOCOM_Com_P=%d | LOCOM_Com_Count=%d | LOCOM_16s=%d | LOCOM_shotgun=%d\n",
  seed,
  length(DA_knockoff), length(DA_Com2seq),
  length(DA_Cauchy_Locom), length(DA_Pool_Locom),
  length(DA_Locom_16s), length(DA_Locom_shotgun)
))

# ==== Export detected taxa names to Excel ==== #
library(writexl)

method_order <- c("DUET-Knockoffs", "Com2seq",
                  "Locom_Com_P", "Locom_Com_Count",
                  "Locom_16s", "Locom_shotgun")

max_len <- max(c(sapply(Taxa_list[method_order], length), 1L))

taxa_df <- as.data.frame(
  lapply(method_order, function(m) {
    taxa <- Taxa_list[[m]]
    c(taxa, rep("", max_len - length(taxa)))
  }),
  stringsAsFactors = FALSE
)
colnames(taxa_df) <- method_order

write_xlsx(taxa_df, path = "Dietary_detected_taxa.xlsx")
cat("Saved selected taxa names to: Dietary_detected_taxa.xlsx\n")
