
<!-- README.md is generated from README.Rmd. Please edit that file -->

# DUET-Knockoffs

## Overview

DUET-Knockoffs is a knockoff-based framework for integrative feature
selection from paired 16S and shotgun microbiome count data. For each
taxon, it tests a union null hypothesis that the taxon is not
conditionally associated with the outcome in at least one platform.
Within each platform, taxon counts are modeled using zero-inflated
negative binomial (ZINB) distributions coupled with a Gaussian copula to
generate knockoff statistics. Platform-specific knockoff contrasts are
then combined via the simultaneous knockoff filter to obtain a single
statistic per taxon, yielding finite-sample false discovery rate (FDR)
control for cross-platform consistent signals.

<p align="center">

<img src="figures/Overview.png" alt="DUET-Knockoffs workflow" width="100%"/>

</p>

## Installation

> **R version (\>= 4.5.1)**

``` r
# install.packages("devtools")
devtools::install_github("dyxstat/DUET-Knockoffs")
```

## Example

An example of using the DUET_Knockoffs function

``` r
library(DUETKnockoffs)

set.seed(2024)

data("W")        # count matrix
data("M")        # library size
data("y")        # control / case group
data("data_x")   # host covariates
data("class_K")  # sequencing platforms

# Run DUET_Knockoffs using differential expression–based test statistic
res.DUET_Knockoffs <- DUET_Knockoffs(
  W = W,
  M = M,
  class_K = class_K,
  y = y,
  data_x = data_x,
  T_var = NULL,          # causal taxa set is NULL
  test_statistic = "DE",
  filter_statistics = 1,
  fdr = 0.1
)

# Detected genera set
genera_idx <- res.DUET_Knockoffs$S
genera_set <- colnames(W)[genera_idx]
genera_set

# [1] "g__Porphyromonas"   "g__Actinomyces"     "g__Neisseria"       "g__Dialister"       "g__Haemophilus"    
#  [6] "g__Capnocytophaga"  "g__Conchiformibius" "g__Prevotella"      "g__Veillonella"     "g__Megasphaera"    
# [11] "g__Campylobacter"   "g__Oribacterium"    "g__Actinobacillus"  "g__Atopobium"       "g__Mogibacterium"  
# [16] "g__Selenomonas"  

filter_stats <- res.DUET_Knockoffs$filter_stat[genera_idx]
names(filter_stats) <- genera_set
filter_stats

# g__Porphyromonas     g__Actinomyces       g__Neisseria       g__Dialister     g__Haemophilus 
#           18.00123           28.99837           38.98954           14.13199           64.86431 
#  g__Capnocytophaga g__Conchiformibius      g__Prevotella     g__Veillonella     g__Megasphaera 
#           15.39400           16.72202          350.72112          182.41365          107.04932 
#   g__Campylobacter    g__Oribacterium  g__Actinobacillus       g__Atopobium   g__Mogibacterium 
#           24.86322           15.11381           24.30750          119.22269           50.84153 
#     g__Selenomonas 
#           32.34956 
```
