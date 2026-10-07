###############################################################################
# 06_classifiers.R
# Comparison with C1A/C1B, CIMP (TCGA pan-genomic analysis) and the
# BUB1B-PINK1 score, including out-of-fold scores (Supplementary Table S6).
###############################################################################

cd <- as.data.frame(colData(acc_data))[, c("barcode", "paper_C1A.C1B", "paper_MethyLevel")]
d <- d %>% dplyr::select(-any_of(c("paper_C1A.C1B", "paper_MethyLevel"))) %>% left_join(cd, by = "barcode")
bub1b <- names(symbol_of)[symbol_of == "BUB1B"][1]; pink1 <- names(symbol_of)[symbol_of == "PINK1"][1]
d$bub1b_pink1 <- (vst_blind[bub1b, ] - vst_blind[pink1, ])[d$barcode]

k <- d[!is.na(d$paper_C1A.C1B) & !is.na(d$paper_MethyLevel), ]
k$C1   <- factor(k$paper_C1A.C1B, levels = c("C1B", "C1A"))
k$CIMP <- factor(k$paper_MethyLevel, levels = c("CIMP-low", "CIMP-intermediate", "CIMP-high"))
Sk <- Surv(k$time, k$event)
cat("Classifier comparison:", nrow(k), "patients,", sum(k$event), "events\n")   # expected 78 / 36
cat("Correlation risk score vs BUB1B-PINK1:", round(cor(k$rz, k$bub1b_pink1), 2), "\n")

cls_rows <- do.call(rbind, lapply(c("C1", "CIMP", "bub1b_pink1"), function(ref) {
  f0 <- coxph(as.formula(paste("Sk ~", ref)), k); f1 <- coxph(as.formula(paste("Sk ~ rz +", ref)), k)
  s <- summary(f1)$conf.int["rz", ]
  data.frame(Reference = ref, C_reference = round(concordance(f0)$concordance, 3),
             C_reference_plus_score = round(concordance(f1)$concordance, 3),
             score_HR_perSD = fmt_ci(s[1], s[3], s[4], 2), LRT_p = signif(anova(f0, f1)[2, 4], 2))
}))
kk <- k[k$C1 == "C1A", ]; s_c1a <- summary(coxph(Surv(time, event) ~ rz, kk))
cat("Within C1A: n =", nrow(kk), "events =", sum(kk$event), "HR/SD =",
    fmt_ci(s_c1a$conf.int[1], s_c1a$conf.int[3], s_c1a$conf.int[4], 2), "\n")

# Out-of-fold scores (5-fold CV, 50 repeats)
set.seed(SEED)
oof <- t(replicate(50, {
  fold <- sample(rep(1:5, length.out = nrow(k))); s <- rep(NA, nrow(k))
  for (ff in 1:5) { tr <- fold != ff
    fit <- coxph(as.formula(paste("Surv(time, event) ~", paste(GENES6, collapse = "+"))), data = k[tr, ])
    s[!tr] <- predict(fit, newdata = k[!tr, ], type = "lp") }
  k$oof <- s; f0 <- coxph(Surv(time, event) ~ C1, k); f1 <- coxph(Surv(time, event) ~ oof + C1, k)
  c(C_oof = concordance(coxph(Surv(time, event) ~ oof, k))$concordance, C_C1 = concordance(f0)$concordance,
    C_both = concordance(f1)$concordance, HR_perSD = unname(exp(coef(f1)["oof"] * sd(s))), LRT_p = anova(f0, f1)[2, 4])
}))
oof_summary <- apply(oof, 2, quantile, c(.5, .025, .975))
print(round(oof_summary, 4))                                          # expected C_oof ~0.749

write.csv(cls_rows, file.path(DIR_RES, "TableS6a_classifiers_TCGA.csv"), row.names = FALSE)
write.csv(as.data.frame(round(oof_summary, 4)), file.path(DIR_RES, "TableS6b_out_of_fold.csv"))
