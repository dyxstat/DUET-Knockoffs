
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DUETknockoff

The DUETknockoff package provides a knockoff-based framework for
integrative feature selection from paired 16S and shotgun microbiome
profiles. DUET-knockoff uses a ZINB–Gaussian copula model to construct
platform-specific knockoffs and, for each taxon, tests a union null that
it is not conditionally associated with the outcome in at least one
platform. Discoveries reflect associations replicated across both
technologies.

## Installation

``` r
# install.packages("devtools")
devtools::install_github("dyxstat/DUET-Knockoffs")
```

## Example

An example of using the DUET_knockoff function

``` r
library(DUETknockoff)

set.seed(42)

data("W")        # count matrix
data("M")        # library size
data("y")        # control / case group
data("data_x")   # host covariates
data("class_K")  # sequencing platforms

# Run DUET_knockoff using differential expression–based test statistic
res.DUETknockoff <- DUET_knockoff(
  W = W,
  M = M,
  class_K = class_K,
  y = y,
  data_x = data_x,
  T_var = NULL,          # causal taxa is NULL
  test_statistic = "DE",
  filter_statistics = 1,
  fdr = 0.1,
)

# Detected genera set
genera_idx <- res.DUETknockoff$S
genera_set <- colnames(W)[genera_idx]
genera_set
```
