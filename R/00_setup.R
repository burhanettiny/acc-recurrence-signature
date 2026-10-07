###############################################################################
# 00_setup.R
# Packages, paths, global settings and helper functions.
# Sourced by run_all.R; all scripts are run in a single R session.
###############################################################################

suppressPackageStartupMessages({
  library(TCGAbiolinks); library(SummarizedExperiment); library(DESeq2)
  library(dplyr); library(tibble); library(tidyr)
  library(clusterProfiler); library(org.Hs.eg.db); library(GO.db)
  library(STRINGdb); library(igraph); library(ggraph)
  library(survival); library(survminer); library(rms); library(cmprsk); library(glmnet)
  library(GEOquery); library(Biobase); library(hgu133plus2.db); library(hgu133a.db)
  library(ggplot2)
})

# ---- Paths ------------------------------------------------------------------
# On Windows, TCGAbiolinks cannot read files from paths that contain non-ASCII
# characters. If the project folder contains such characters, set GDC_DIR to an
# ASCII path before running, e.g. Sys.setenv(GDC_DIR = "C:/GDCdata").
DIR_BIOTAB <- file.path("data", "raw_clinical_biotab")
DIR_RES    <- "results"
DIR_FIG    <- "figures"
DIR_MOD    <- "models"
GDC_DIR    <- Sys.getenv("GDC_DIR", unset = file.path("data", "GDCdata"))
for (d in c(DIR_RES, DIR_FIG, DIR_MOD, GDC_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

SEED <- 2026
options(timeout = 1200)   # large downloads (GDC, STRING, GEO)
DIR_STRING <- Sys.getenv("STRING_DIR", unset = file.path("data", "STRINGdb_cache"))
dir.create(DIR_STRING, showWarnings = FALSE, recursive = TRUE)
# The ridge coefficients reported in the manuscript are archived in models/.
# With USE_ARCHIVED_MODELS = TRUE the archived coefficients are used for all
# downstream analyses (the models are also refitted and compared, see 04_).
USE_ARCHIVED_MODELS <- TRUE

GENES7 <- c("NDC80", "BUB1", "TYMS", "ASPM", "NCAPH", "CDK6", "SPAG5")
GENES6 <- setdiff(GENES7, "CDK6")
GENES5 <- c("NDC80", "BUB1", "ASPM", "NCAPH", "SPAG5")   # GO:0051301 rule without CDK6

# ---- Helpers ----------------------------------------------------------------
num <- function(x) suppressWarnings(as.numeric(x))

# Read one BCR Biotab clinical table; the first two data rows are CDE descriptors.
read_biotab <- function(table_name) {
  f <- list.files(DIR_BIOTAB, full.names = TRUE)
  f <- f[endsWith(f, paste0("_clinical_", table_name, ".txt"))]
  stopifnot(length(f) == 1)
  x <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE, quote = "", na.strings = "")
  names(x) <- make.unique(names(x))
  x[-c(1, 2), , drop = FALSE]
}

# Risk score from within-cohort z-scored expression and per-SD coefficients b
score_z <- function(M, b) as.vector(scale(as.matrix(M[, names(b), drop = FALSE])) %*% b)

fmt_ci <- function(est, lo, hi, d = 3) sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), est, lo, hi)

# One-row summary of a continuous score: C-index (95% CI), HR per SD (95% CI), p, PH p
cox_row <- function(time, event, score, label, cohort = NA) {
  y <- Surv(time, event); f <- coxph(y ~ scale(score))
  cc <- concordance(f); se <- sqrt(cc$var); ci <- summary(f)$conf.int
  data.frame(Model = label, Cohort = cohort, n = length(time), events = sum(event),
             C = fmt_ci(cc$concordance, cc$concordance - 1.96 * se, cc$concordance + 1.96 * se),
             HR_perSD = fmt_ci(ci[1], ci[3], ci[4], 2),
             p = signif(summary(f)$coefficients[5], 3),
             PH_p = signif(cox.zph(f)$table[1, 3], 2))
}

# Mean of all probe sets per gene (rule used for all microarray cohorts)
gene_matrix <- function(ex, annot_db, genes, samples) {
  pm <- AnnotationDbi::select(annot_db, keys = genes, keytype = "SYMBOL", columns = "PROBEID")
  pm <- pm[pm$PROBEID %in% rownames(ex), ]
  sapply(genes, function(g) colMeans(ex[pm$PROBEID[pm$SYMBOL == g], samples, drop = FALSE]))
}

set.seed(SEED)
cat("Setup complete. R", R.version$major, ".", R.version$minor, "\n")
