## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>"
)

## ----setup--------------------------------------------------------------------
library(SNPkit)
library(snpStats)
library(methods)

## ----build--------------------------------------------------------------------
set.seed(123)

raw_mat <- matrix(
  as.raw(sample(1:3, 100, replace = TRUE)),
  nrow = 10, ncol = 10
)
rownames(raw_mat) <- paste0("ind", 1:10)
colnames(raw_mat) <- paste0("snp", 1:10)

geno <- new("SnpMatrix", raw_mat)

map <- data.frame(
  Name       = colnames(geno),
  Chromosome = rep(1, 10),
  Position   = seq_len(10),
  stringsAsFactors = FALSE
)

snp_data <- new(
  "SNPDataLong",
  geno      = geno,
  map       = map,
  path      = tempfile(),
  xref_path = "chip1"
)

snp_data

## ----summary------------------------------------------------------------------
s <- summary(snp_data)
s$n_individuals
s$n_snps
s$prop_missing
print(s)

## ----qc-----------------------------------------------------------------------
filtered <- qcSNPs(
  snp_data,
  min_snp_cr   = 0.8,
  min_maf      = 0.05,
  snp_mono     = TRUE,
  no_position  = TRUE,
  action       = "filter"
)
filtered

## ----export-------------------------------------------------------------------
out_dir <- file.path(tempdir(), "snpkit_demo")
dir.create(out_dir, showWarnings = FALSE)

savePlink(
  filtered,
  path       = out_dir,
  name       = "demo",
  run_plink  = FALSE,
  chunk_size = 5
)
list.files(out_dir, pattern = "demo")

