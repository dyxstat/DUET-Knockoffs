
# Expected directory layout (relative to this script):
# - 1_ORIGINS_processing.R  (creates meta_common, otu_tab_commonSam, shotgun_tab_commonSam, etc.)
# - 2_ORIGINS_analOverlappingData_Figure2-6_Table1-2.R (not required to run but kept for reference)
# - 3_ORIGINS_analFullData_Figure7_Table1-3.R (not required to run but kept for reference)
# - ORIGINS raw data in ./ID11808/... as referenced by 1_ORIGINS_processing.R

# Dependencies:
# - R >= 4.1
# - Packages: DUETknockoff, knockoff, dirmult, readxl, biomformat, GUniFrac, vegan, readr, 
#          LOCOM, matrixStats, prodlim, parallel
#
# install.packages("devtools")
devtools::install_github("dyxstat/DUET-Knockoffs")

suppressPackageStartupMessages({
  library(DUETknockoff)       
  library(dirmult)
  library(readxl)
  library(readr)
  library(biomformat)
  library(GUniFrac)
  library(vegan)
  library(matrixStats)
  library(LOCOM)
})

set.seed(42)

# ----------------------------------------------------------------------------
# 0) Load & prepare ORIGINS objects using the user's processing script
# ----------------------------------------------------------------------------
message("Sourcing ORIGINS_processing.R to load ORIGINS data...")
source("ORIGINS_processing.R")

# Required objects created by 1_ORIGINS_processing.R:
# - meta_common               : metadata for overlapping 16S/SMS samples (n x p)
# - otu_tab_commonSam         : 16S counts (n x G1) for overlapping samples
# - shotgun_tab_commonSam     : SMS counts (n x G2) for overlapping samples
# - filter_thresh             : prevalence threshold used in upstream scripts (default 0.2)
stopifnot(exists("meta_common"), exists("otu_tab_commonSam"), exists("shotgun_tab_commonSam"))

results_dir <- "results"
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

# ----------------------------------------------------------------------------
# 1) Align taxa across 16S and SMS: restrict to common genera, then filter by prevalence
# ----------------------------------------------------------------------------
common_taxa <- intersect(colnames(otu_tab_commonSam), colnames(shotgun_tab_commonSam))
if (length(common_taxa) < 5) {
  stop("Fewer than 5 common taxa between 16S and SMS; check preprocessing or taxon names.")
}

otu_common     <- otu_tab_commonSam[, common_taxa, drop = FALSE]
shotgun_common <- shotgun_tab_commonSam[, common_taxa, drop = FALSE]

pooled_counts <- otu_common + shotgun_common
prev <- colMeans(pooled_counts > 0)
keep <- which(prev >= if (exists("filter_thresh")) filter_thresh else 0.2)
if (length(keep) < 5) stop("Too few taxa after prevalence filter; try lowering filter_thresh.")
otu_common     <- otu_common[, keep, drop = FALSE]
shotgun_common <- shotgun_common[, keep, drop = FALSE]

# ----------------------------------------------------------------------------
# 2) Build covariates X using meta_common.
#    We'll use age (numeric), BMI (numeric), sex, race/ethnicity, smoking.
#    Factors are one-hot encoded via model.matrix to be robust.
# ----------------------------------------------------------------------------
covars <- c("host_age","sex","bmi","cigcurr")

# normalize string-missing tokens just for these four columns
na_tok <- c("not provided","missing","","na","n/a")
meta_common[covars] <- lapply(meta_common[covars], function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    z <- tolower(trimws(x))
    x[z %in% na_tok] <- NA
  }
  x
})

idx <- complete.cases(meta_common[, covars, drop = FALSE])

meta_common    <- meta_common[idx, , drop = FALSE]
otu_common     <- otu_common[idx, , drop = FALSE]
shotgun_common <- shotgun_common[idx, , drop = FALSE]

# build design matrix: scale numeric only
X <- cbind(
  host_age = as.numeric(scale(as.numeric(meta_common$host_age))),
  sex      = as.numeric(meta_common$sex),
  bmi      = as.numeric(scale(as.numeric(meta_common$bmi))),
  cigcurr  = as.numeric(meta_common$cigcurr)
)

data_x_1 <- as.data.frame(X)  # covariates for 16S block
data_x_2 <- as.data.frame(X)  # covariates for SMS block

# >>> Combine covariates for both blocks in the SAME row order as W
data_x <- rbind(data_x_1, data_x_2)


# ----------------------------------------------------------------------------
# 3) Build DUETknockoff inputs
# ----------------------------------------------------------------------------
n_sam <- nrow(otu_common)
stopifnot(n_sam == nrow(shotgun_common))

M1 <- rowSums(otu_common)     # library sizes
M2 <- rowSums(shotgun_common)

# Combine rows: 16S on top of SMS (as in the example)
W1 <- as.data.frame(otu_common)
W2 <- as.data.frame(shotgun_common)
W  <- rbind(W1, W2)

# Combine covariates and library sizes in the same row order
data_x <- rbind(data_x_1, data_x_2)
M      <- c(M1, M2)

# Data-source labels
class_K <- factor(c(rep(1, n_sam), rep(2, n_sam)))

# Binary outcome (1 = Not; 2 = Prediabetes) and then stacked twice (one per data source)
y_single <- ifelse(meta_common$prediabetes %in% c("1","Prediabetes"), 2, 1)
y <- c(y_single, y_single)

# ----------------------------------------------------------------------------
# 4) Run DUET_knockoff
# ----------------------------------------------------------------------------
fdr_target <- 0.2  # adjust as needed

write.csv(W, file.path(results_dir, "W.csv"))
write.csv(class_K, file.path(results_dir, "class_K.csv"))
write.csv(data_x, file.path(results_dir, "data_x.csv"))
write.csv(M, file.path(results_dir, "M.csv"))
write.csv(y, file.path(results_dir, "y.csv"))
write.csv(metaX, file.path(results_dir, "metaX.csv"))


# T_var is unknown for real data; omit it. (It is only for simulation-based power.)
T_var = c(1,2,3)
library(DUETknockoff)
library(irlba)
fit <- DUET_knockoff(
  W = W,
  class_K = class_K,
  data_x  = data_x,
  M       = M,
  y       = y,
  T_var = T_var,
  fdr     = fdr_target,
  test_statistic = "DE",
  filter_statistics = 3,
  offset  = 1
)

# ----------------------------------------------------------------------------
# 5) Save outputs
# ----------------------------------------------------------------------------
saveRDS(fit, file = file.path(results_dir, "DUETknockoff_full_result.RDS"))

# Try to extract selected taxa if 'S' exists (indices of discoveries)
sel_taxa <- NULL
if (!is.null(fit$S)) {
  # fit$S may be indices relative to columns of W (common taxa). Deduplicate if needed.
  sel_idx <- sort(unique(as.integer(fit$S)))
  sel_idx <- sel_idx[sel_idx >= 1 & sel_idx <= ncol(W)]
  sel_taxa <- colnames(W)[sel_idx]
}

# Write a CSV with all taxa and a "selected" flag (if available)
all_taxa <- colnames(W)
df_out <- data.frame(
  taxon = all_taxa,
  selected = if (is.null(sel_taxa)) FALSE else all_taxa %in% sel_taxa,
  stringsAsFactors = FALSE
)
write.csv(df_out, file.path(results_dir, "DUETknockoff_selected_taxa.csv"), row.names = FALSE)

# Summary
sink(file.path(results_dir, "DUETknockoff_summary.txt"))
cat("[DUETknockoff] ORIGINS overlapping-sample run\n")
cat("Date: ", as.character(Sys.time()), "\n", sep = "")
cat("FDR target: ", fdr_target, "\n", sep = "")
cat("Samples (per source): ", n_sam, "\n", sep = "")
cat("Common taxa used: ", ncol(W), "\n", sep = "")
if (!is.null(sel_taxa)) {
  cat("Selected taxa (count): ", length(sel_taxa), "\n", sep = "")
  cat("First few: ", paste(utils::head(sel_taxa, 20), collapse = ", "), "\n", sep = "")
} else {
  cat("Selected taxa: not provided in fit (fit$S missing)\n")
}
if (!is.null(fit$FDRPower)) {
  cat("FDRPower (if present):\n")
  print(fit$FDRPower)
}
sink()

message("[DUETknockoff] Done. See ./results/ for outputs.")

