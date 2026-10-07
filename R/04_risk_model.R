###############################################################################
# 04_risk_model.R
# Hub-gene survival analysis, unpenalized and ridge-penalized Cox models,
# bootstrap optimism, nested cross-validation (Figure 3, Supplementary Table S7).
###############################################################################

# ---- Analysis dataset --------------------------------------------------------
gene_id_map <- sig_genes_annotated %>% filter(SYMBOL %in% GENES7) %>% distinct(SYMBOL, gene_id)
expr_sel <- as.data.frame(t(vst_matrix[gene_id_map$gene_id, , drop = FALSE]))
colnames(expr_sel) <- gene_id_map$SYMBOL[match(colnames(expr_sel), gene_id_map$gene_id)]
expr_sel$barcode <- rownames(expr_sel)
expr_sel <- expr_sel %>% left_join(sample_info[, c("barcode", "submitter_id")], by = "barcode")
d <- survival_df %>% inner_join(expr_sel, by = "submitter_id") %>%
  left_join(clin, by = "submitter_id") %>% arrange(barcode)
cat("Analysis dataset:", nrow(d), "patients,", sum(d$event), "events\n")   # expected 79 / 37
S <- Surv(d$time, d$event)

# ---- Figure 3: Kaplan-Meier curves per gene (median split) -------------------
km_list <- lapply(GENES6, function(gn) {
  df <- data.frame(t = d$time, e = d$event,
                   grp = factor(ifelse(d[[gn]] >= median(d[[gn]]), "High", "Low"), levels = c("Low", "High")))
  ggsurvplot(survfit(Surv(t, e) ~ grp, data = df), data = df, pval = TRUE, risk.table = FALSE,
             title = gn, xlab = "Days", ylab = "Recurrence-free survival", legend.title = "",
             legend.labs = c("Low", "High"), palette = c("#2166AC", "#B2182B"))
})
png(file.path(DIR_FIG, "Figure3_KM_genes.png"), width = 4200, height = 2600, res = 300)
print(arrange_ggsurvplots(km_list, ncol = 3, nrow = 2, print = FALSE)); dev.off()

# ---- Univariate and unpenalized multivariable Cox ----------------------------
uni <- do.call(rbind, lapply(GENES7, function(gn) {
  s <- summary(coxph(as.formula(paste("S ~", gn)), data = d))
  data.frame(gene = gn, HR = s$conf.int[1], lower = s$conf.int[3], upper = s$conf.int[4], p = s$coefficients[5])
}))
print(uni); write.csv(uni, file.path(DIR_RES, "Cox_univariate_7genes.csv"), row.names = FALSE)

m6 <- coxph(as.formula(paste("S ~", paste(GENES6, collapse = "+"))), data = d)
m7 <- coxph(as.formula(paste("S ~", paste(GENES7, collapse = "+"))), data = d)
coefs_6gene <- coef(m6); coefs_7gene <- coef(m7)

# ---- Bootstrap optimism with coefficient re-estimation (B = 500) ------------
set.seed(SEED)
dd6 <- datadist(d[, c(GENES6, "time", "event")]); options(datadist = "dd6")
f6  <- cph(as.formula(paste("Surv(time, event) ~", paste(GENES6, collapse = "+"))), data = d, x = TRUE, y = TRUE)
val6 <- validate(f6, B = 500)
boot_summary <- c(C_apparent = val6["Dxy", "index.orig"] / 2 + 0.5,
                  C_corrected = val6["Dxy", "index.corrected"] / 2 + 0.5,
                  slope = val6["Slope", "index.corrected"])
print(round(boot_summary, 3))                                          # expected 0.790 / 0.756 / 0.82

# ---- Ridge-penalized Cox models -----------------------------------------------
fit_ridge <- function(genes) {
  cv <- cv.glmnet(scale(as.matrix(d[, genes])), S, family = "cox", alpha = 0, nfolds = 10)
  setNames(as.vector(coef(cv, s = "lambda.min")), genes)
}
b6_refit <- fit_ridge(GENES6)                       # continues the RNG stream after validate()

# ---- Nested cross-validation of the full pipeline (5-fold x 20) --------------
V <- vst_blind[, d$barcode]; sym_V <- symbol_of[rownames(V)]
cvC <- c(); picked <- list()
for (r in 1:20) {
  fold <- sample(rep(1:5, length.out = ncol(V)))
  for (ff in 1:5) {
    tr <- fold != ff; te <- !tr; Vt <- V[, tr]; g <- d$event[tr]
    lfc  <- rowMeans(Vt[, g == 1]) - rowMeans(Vt[, g == 0])
    se   <- sqrt(apply(Vt[, g == 1], 1, var) / sum(g == 1) + apply(Vt[, g == 0], 1, var) / sum(g == 0))
    pv   <- 2 * pt(-abs(lfc / se), df = sum(tr) - 2); padj <- p.adjust(pv, "BH")
    cand <- which(padj < 0.05 & abs(lfc) > 1); cand <- cand[order(pv[cand])][1:min(300, length(cand))]
    ytr  <- Surv(d$time[tr], d$event[tr])
    z    <- sapply(cand, function(i) summary(coxph(ytr ~ Vt[i, ]))$coefficients[4])
    pick <- cand[order(-abs(z))][1:6]; picked[[length(picked) + 1]] <- sym_V[pick]
    fit  <- cv.glmnet(t(Vt[pick, ]), ytr, family = "cox", alpha = 0, nfolds = 5)
    lp   <- as.vector(t(V[pick, te]) %*% as.vector(coef(fit, s = "lambda.min")))
    cvC  <- c(cvC, concordance(Surv(d$time[te], d$event[te]) ~ lp, reverse = TRUE)$concordance)
  }
  cat("  nested CV repeat", r, "of 20 done\n")
}
nested_summary <- c(mean = mean(cvC), median = median(cvC), quantile(cvC, c(.025, .975)))
print(round(nested_summary, 3))                                        # expected mean ~0.70
print(sort(table(unlist(picked)), decreasing = TRUE)[1:10])

set.seed(SEED); b7_refit <- fit_ridge(GENES7)
set.seed(SEED); b5_refit <- fit_ridge(GENES5)

# ---- Coefficients used downstream ---------------------------------------------
arch_file <- file.path(DIR_MOD, "archived_ridge_coefficients.rds")
if (USE_ARCHIVED_MODELS && file.exists(arch_file)) {
  arch <- readRDS(arch_file)
  b6 <- arch$b6; b7 <- arch$b7; b5 <- arch$b5
  cat("Max |difference| refit vs archived (b6, b7, b5):",
      signif(c(max(abs(b6 - b6_refit)), max(abs(b7 - b7_refit)), max(abs(b5 - b5_refit))), 3), "\n")
} else {
  b6 <- b6_refit; b7 <- b7_refit; b5 <- b5_refit
}
print(round(b6, 3))   # NDC80 0.043, BUB1 0.255, TYMS 0.363, ASPM 0.241, NCAPH 0.152, SPAG5 0.138

d$ridge_score <- score_z(d, b6)
d$rz <- as.vector(scale(d$ridge_score))

tabS7 <- data.frame(gene = GENES6, beta_unpenalized = coef(m6), HR = summary(m6)$conf.int[, 1],
                    lower = summary(m6)$conf.int[, 3], upper = summary(m6)$conf.int[, 4],
                    p = summary(m6)$coefficients[, 5], beta_ridge_perSD = b6[GENES6])
write.csv(tabS7, file.path(DIR_RES, "TableS7_coefficients.csv"), row.names = FALSE)
saveRDS(list(b6 = b6, b7 = b7, b5 = b5, coefs_6gene = coefs_6gene, coefs_7gene = coefs_7gene,
             boot = boot_summary, nested = nested_summary, nested_C = cvC),
        file.path(DIR_MOD, "models_this_run.rds"))
