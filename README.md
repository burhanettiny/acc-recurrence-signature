# Recurrence-associated cell-cycle gene signature in adrenocortical carcinoma

Code to reproduce all analyses of:

> Yalçınkaya B, *et al.* Identification of a Recurrence-Associated Cell-Cycle Gene Signature in Adrenocortical Carcinoma: Development in TCGA and External Evaluation in Adult and Pediatric Cohorts. *[Journal]*, [year]. doi: [to be added]

A six-gene signature (*NDC80, BUB1, TYMS, ASPM, NCAPH, SPAG5*) associated with recurrence-free survival after surgical resection was derived in TCGA-ACC with ridge-penalized Cox regression and evaluated in three GEO cohorts (GSE10927, GSE19750, GSE76021).

## Repository structure

```
run_all.R                    runs the complete analysis in one R session
R/00_setup.R                 packages, paths, helper functions
R/01_tcga_data.R             TCGA-ACC RNA-seq, BCR Biotab clinical data, recurrence definition, RFS
R/02_deseq2_enrichment.R     DESeq2, stage-adjusted DESeq2, GO/KEGG (Figure 1, Table S1)
R/03_ppi_hub_selection.R     STRING network, centrality, hub-gene selection (Figure 2, Table S2)
R/04_risk_model.R            Cox models, ridge, bootstrap, nested cross-validation (Figure 3, Table S7)
R/05_clinical_nomogram.R     Table 1, multivariable and restricted analyses, competing risks,
                             nomogram and calibration (Tables 1-2, S3; Figures 4-5)
R/06_classifiers.R           C1A/C1B, CIMP and BUB1B-PINK1 comparison (Table S6)
R/07_external_validation.R   GEO cohorts (Table 3, Tables S4/S5/S8/S9, Figure 6)
data/raw_clinical_biotab/    TCGA BCR Biotab clinical files used in the study (+ GDC MANIFEST)
data/tcga_followup_snapshot.csv   follow-up data used in the study (see below)
models/archived_ridge_coefficients.rds   ridge coefficients reported in the manuscript
tools/export_from_workspace.R    one-time export used by the authors (not needed to reproduce)
```

Outputs are written to `results/` (tables, `sessionInfo.txt`) and `figures/`.

## Requirements

R ≥ 4.5 with Bioconductor 3.22. The analyses were run with R 4.5.2, DESeq2 1.50.2, TCGAbiolinks 2.38.0, clusterProfiler 4.18.4, STRINGdb 2.22.0, igraph 2.3.2, survival 3.8.6, survminer 0.5.2, rms 8.1.1, cmprsk 2.2.12, glmnet 5.0, GEOquery 2.78.0, msigdbr 26.1.1, hgu133plus2.db 3.13.0 and hgu133a.db 3.13.0 (full list in `results/sessionInfo.txt`).

```r
install.packages(c("BiocManager", "dplyr", "tibble", "tidyr", "survival", "survminer",
                   "rms", "cmprsk", "glmnet", "igraph", "ggraph", "ggplot2", "polspline"))
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", "clusterProfiler",
                       "org.Hs.eg.db", "GO.db", "STRINGdb", "GEOquery", "Biobase",
                       "hgu133plus2.db", "hgu133a.db", "EnhancedVolcano"))
```

## Running

From the repository root:

```bash
Rscript run_all.R
```

The first run downloads TCGA-ACC RNA-seq data (~1 GB) from the GDC and three GEO series; an internet connection is required. On Windows, TCGAbiolinks cannot read files from paths containing non-ASCII characters; if necessary, set an ASCII download folder before running, e.g. `Sys.setenv(GDC_DIR = "C:/GDCdata")`.

## Reproducibility notes

* **Clinical data.** GDC clinical records are updated over time. The BCR Biotab files and the follow-up values used in the study are bundled in `data/` and are used in preference to newly downloaded records.
* **Ridge coefficients.** The coefficients reported in the manuscript are archived in `models/`. The models are also refitted in `04_risk_model.R`, and the maximum difference between refitted and archived coefficients is printed. Set `USE_ARCHIVED_MODELS <- FALSE` in `00_setup.R` to use the refitted coefficients instead.
* **Resampling.** Bootstrap and cross-validation estimates (optimism-corrected C-index, nested cross-validation, calibration, confidence intervals of ΔC) depend on the random-number stream and may differ slightly (typically in the third decimal) between platforms and package versions.
* **Expected key values** are printed as comments next to the corresponding commands (e.g. 79 patients / 37 events; 2,101 DEGs; ridge C-index 0.773).

## Data sources

* TCGA-ACC: Genomic Data Commons, https://portal.gdc.cancer.gov (project TCGA-ACC).
* GEO: GSE10927, GSE19750, GSE76021, https://www.ncbi.nlm.nih.gov/geo.

## Citation

If you use this code, please cite the article above and the archived release: [https://doi.org/10.5281/zenodo.23213384].

## License

MIT (see `LICENSE`).
