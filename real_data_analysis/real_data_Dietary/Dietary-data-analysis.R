

# Run DUETknockoff on the second real dataset ("Diet")

### ─────────────────────────────────────────────────────────────────────────────
### 1) Source DIET processing script (creates the objects)
### ─────────────────────────────────────────────────────────────────────────────
source("Diet_processing.R")

if (!exists("meta")) {
  for (nm in c("meta_diet","meta_common","meta")) if (exists(nm)) { meta <- get(nm); break }
}
if (!exists("otu")) {
  for (nm in c("otu_tab_commonSam","otu_common","otu_diet","otu_tab_diet_commonSam")) if (exists(nm)) { otu <- get(nm); break }
}
if (!exists("sms")) {
  for (nm in c("shotgun_tab_commonSam","shotgun_common","sms_diet","shotgun_tab_diet_commonSam")) if (exists(nm)) { sms <- get(nm); break }
}
stopifnot(exists("meta"), exists("otu"), exists("sms"))
cat("Loaded objects:\n"); print(list(meta=dim(meta), otu=dim(otu), sms=dim(sms)))

### ─────────────────────────────────────────────────────────────────────────────
### 2) Align taxa across 16S and SMS: restrict to common genera, then filter by prevalence
### ─────────────────────────────────────────────────────────────────────────────
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



### ─────────────────────────────────────────────────────────────────────────────
### 2) Choose covariates and drop rows with missing values 
### ─────────────────────────────────────────────────────────────────────────────
na_tok <- c("not provided","missing","","na","n/a","nan")
v <- meta_common$phylogency
v <- if (is.factor(v)) as.character(v) else v
v[tolower(trimws(v)) %in% na_tok] <- NA
meta_common$phylogency <- v
idx <- complete.cases(meta_common[, "phylogency", drop = FALSE])
meta_common    <- meta_common[idx, , drop = FALSE]
otu_common     <- otu_common[idx, , drop = FALSE]
shotgun_common <- shotgun_common[idx, , drop = FALSE]


## One-hot encode phylogency (no intercept → 1 column per level)
meta_common$phylogency <- factor(meta_common$phylogency)
X <- model.matrix(~ phylogency - 1, data = meta_common)
colnames(X) <- sub("^phylogency", "", colnames(X))


data_x_1 <- as.data.frame(X)
data_x_2 <- as.data.frame(X)


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

# Binary outcome (1 = Not; 2 = Folivore) and then stacked twice (one per data source)
y_single <- ifelse(meta_common$diet %in% c("1","Folivore"), 2, 1)
y <- c(y_single, y_single)


# ----------------------------------------------------------------------------
# 4) Run ZINB-SK
# ----------------------------------------------------------------------------
fdr_target <- 0.2  # adjust as needed

write.csv(W, "W.csv")
write.csv(class_K, "class_K.csv")
write.csv(data_x, "data_x.csv")
write.csv(M, "M.csv")
write.csv(y, "y.csv")
write.csv(metaX, "metaX.csv")


# T_var is unknown for real data; omit it. (It is only for simulation-based power.)
T_var = c(1,2,3)
library(DUETknockoff)
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
print(tm); names(fit)

# Extract selected taxa if available
sel_idx <- NULL
if (!is.null(fit$S)) {
  sel_idx <- sort(unique(as.integer(fit$S)))
  sel_idx <- sel_idx[sel_idx >= 1 & sel_idx <= ncol(W)]
}
sel_taxa <- if (length(sel_idx)) colnames(W)[sel_idx] else character(0)
cat("Selected taxa:", length(sel_taxa), "\n")
if (length(sel_taxa)) print(utils::head(sel_taxa, 20))

### ─────────────────────────────────────────────────────────────────────────────
### 5) Save outputs
### ─────────────────────────────────────────────────────────────────────────────
outdir <- "results_diet"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

saveRDS(fit, file = file.path(outdir, "DUETknockoff_full_result_diet.RDS"))
write.csv(
  data.frame(taxon = colnames(W), selected = colnames(W) %in% sel_taxa, stringsAsFactors = FALSE),
  file = file.path(outdir, "DUETknockoff_selected_taxa_diet.csv"),
  row.names = FALSE
)
sink(file.path(outdir, "DUETknockoff_summary_diet.txt"))
cat("Date: ", as.character(Sys.time()), "\n", sep = "")
cat("Outcome column: ", outcome_col, "\n", sep = "")
cat("Covariates: ", paste(covars, collapse=", "), "\n", sep = "")
cat("FDR target: ", fdr_target, "\n", sep = "")
cat("Samples per source: ", n, "\n", sep = "")
cat("Taxa used: ", ncol(W), "\n", sep = "")
cat("Selected taxa: ", length(sel_taxa), "\n", sep = "")
if (length(sel_taxa)) cat("First few: ", paste(utils::head(sel_taxa, 20), collapse = ", "), "\n", sep = "")
sink()
cat("Saved outputs to:", normalizePath(outdir), "\n")
