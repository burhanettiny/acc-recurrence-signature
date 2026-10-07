###############################################################################
# tools/export_from_workspace.R  (run ONCE by the authors, not needed by users)
# Exports, from the original analysis session, the two items that cannot be
# regenerated identically from public sources:
#   1. the follow-up snapshot (GDC clinical records change over time);
#   2. the ridge coefficients reported in the manuscript.
# Run from the folder that contains the original revision files.
###############################################################################
e <- new.env()
load("ACC_analysis_FULL_WORKSPACE.RData", envir = e)
snap <- e$surv_expr_df[, c("submitter_id", "vital_status", "days_to_last_follow_up", "days_to_death")]
write.csv(snap, "tcga_followup_snapshot.csv", row.names = FALSE)

ri <- readRDS("ridge_model.rds")          # six-gene ridge model (b, lambda, center, scale)
s8 <- readRDS("step8.rds")                # contains b7 (seven-gene ridge)
d  <- readRDS("d_revision_final.rds")     # analysis dataset used for the revision
library(glmnet); library(survival)
g5 <- c("NDC80", "BUB1", "ASPM", "NCAPH", "SPAG5")
set.seed(2026)
b5 <- setNames(as.vector(coef(cv.glmnet(scale(as.matrix(d[, g5])), Surv(d$time, d$event),
                                         family = "cox", alpha = 0, nfolds = 10), s = "lambda.min")), g5)
saveRDS(list(b6 = ri$b, b7 = s8$b7, b5 = b5, lambda6 = ri$lambda), "archived_ridge_coefficients.rds")
print(lapply(list(b6 = ri$b, b7 = s8$b7, b5 = b5), round, 4))
