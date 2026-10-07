###############################################################################
# run_all.R  -  reproduces all analyses of the manuscript
# Usage (from the repository root):  Rscript run_all.R
###############################################################################
t0 <- Sys.time()
for (f in c("00_setup.R", "01_tcga_data.R", "02_deseq2_enrichment.R", "03_ppi_hub_selection.R",
            "04_risk_model.R", "05_clinical_nomogram.R", "06_classifiers.R", "07_external_validation.R")) {
  cat("\n==========", f, "==========\n")
  source(file.path("R", f), echo = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path("results", "sessionInfo.txt"))
cat("\nFinished in", round(difftime(Sys.time(), t0, units = "mins"), 1), "minutes\n")
