###############################################################################
# 05_clinical_nomogram.R
# Table 1, multivariable models and added value, restricted analyses (Table S3),
# competing risks, proportional hazards, nomogram and calibration (Figures 4-5),
# Table 2.
###############################################################################

d$grp <- factor(ifelse(d$event == 1, "Recurrent", "Non-recurrent"), levels = c("Non-recurrent", "Recurrent"))

# ---- Table 1 -----------------------------------------------------------------
cat_rows <- function(v, lab) {
  t <- table(d[[v]], d$grp); p <- fisher.test(t)$p.value
  out <- data.frame(Variable = paste0(lab, ": ", rownames(t)),
                    `Non-recurrent` = sprintf("%d (%.1f%%)", t[, 1], 100 * t[, 1] / sum(t[, 1])),
                    Recurrent = sprintf("%d (%.1f%%)", t[, 2], 100 * t[, 2] / sum(t[, 2])),
                    p = c(signif(p, 3), rep(NA, nrow(t) - 1)), check.names = FALSE)
  miss <- tapply(is.na(d[[v]]), d$grp, sum)
  rbind(out, data.frame(Variable = paste0(lab, ": missing"), `Non-recurrent` = miss[1], Recurrent = miss[2],
                        p = NA, check.names = FALSE))
}
num_row <- function(v, lab) {
  q <- tapply(d[[v]], d$grp, function(x) sprintf("%.0f [%.0f-%.0f]", median(x, na.rm = TRUE),
                                                 quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE)))
  data.frame(Variable = lab, `Non-recurrent` = q[1], Recurrent = q[2],
             p = signif(wilcox.test(d[[v]] ~ d$grp, exact = FALSE)$p.value, 3), check.names = FALSE)
}
tab1 <- rbind(num_row("age", "Age, years"), cat_rows("sex", "Sex"), cat_rows("stageS", "ENSAT stage"),
              cat_rows("hormone", "Hormone status"), cat_rows("residual_tumor", "Resection"),
              cat_rows("mitotane_adjuvant", "Adjuvant mitotane"), num_row("time", "Time to event or censoring, days"))
rownames(tab1) <- NULL
write.csv(tab1, file.path(DIR_RES, "Table1_baseline.csv"), row.names = FALSE)
cat("First NTE type among recurrent patients:\n"); print(table(d$first_nte_type[d$event == 1], useNA = "ifany"))

# ---- Multivariable models and added value ------------------------------------
ds <- d[!is.na(d$stageS), ]
f_stage <- coxph(Surv(time, event) ~ rz + stageS, data = ds)
dm <- d[complete.cases(d[, c("rz", "stageS", "age", "hormone", "nonR0")]), ]
full <- coxph(Surv(time, event) ~ rz + stageS + age + hormone + nonR0, data = dm)
clinm <- coxph(Surv(time, event) ~ stageS + age + hormone + nonR0, data = dm)
set.seed(SEED)
dC <- replicate(500, { b <- dm[sample(nrow(dm), replace = TRUE), ]
  concordance(coxph(Surv(time, event) ~ rz + stageS + age + hormone + nonR0, data = b))$concordance -
  concordance(coxph(Surv(time, event) ~ stageS + age + hormone + nonR0, data = b))$concordance })
added <- c(C_clin = concordance(clinm)$concordance, C_full = concordance(full)$concordance,
           dC = concordance(full)$concordance - concordance(clinm)$concordance,
           dC_lo = unname(quantile(dC, .025)), dC_hi = unname(quantile(dC, .975)),
           LRT_p = anova(clinm, full)[2, 4])
print(summary(full)$conf.int); print(round(added, 4))       # expected HR/SD 3.15; C 0.750 -> 0.812

# ---- Restricted / sensitivity analyses (Table S3) ---------------------------
r_row <- function(dd, lab) cox_row(dd$time, dd$event, dd$rz, lab)[, c("Model", "n", "events", "C", "HR_perSD", "p")]
tabS3 <- rbind(
  r_row(d, "All patients"),
  r_row(d[d$residual_tumor %in% "R0", ], "R0 resection only"),
  r_row(d[!grepl("IV", d$stage_raw), ], "Stage IV excluded"),
  r_row(d[d$residual_tumor %in% "R0" & !grepl("IV", d$stage_raw), ], "R0 and stage I-III"),
  r_row(d[d$stageS %in% "Early (I-II)", ], "Stage I-II only"),
  r_row(d[!(d$event == 0 & d$time < 730), ], "Non-recurrent follow-up >= 2 years"),
  r_row(d[!(d$event == 0 & d$time < 1095), ], "Non-recurrent follow-up >= 3 years"),
  r_row(d[!(d$first_nte_type %in% "New Primary Tumor"), ], "New primary tumour excluded"),
  r_row(d[!is.na(d$age) & d$age >= 18, ], "Adults only (age >= 18)"))
f_res <- coxph(Surv(time, event) ~ rz + nonR0, data = d)
print(tabS3); print(summary(f_res)$conf.int)
write.csv(tabS3, file.path(DIR_RES, "TableS3_sensitivity.csv"), row.names = FALSE)

# ---- Competing risks and proportional hazards ------------------------------
status <- ifelse(d$event == 1, 1, ifelse(d$vital_status == "Dead", 2, 0))
fg <- crr(d$time, status, cov1 = matrix(d$rz), failcode = 1, cencode = 0)
cat("Competing deaths:", sum(status == 2), "\n"); print(summary(fg)$conf.int)   # expected sHR 3.05
f_rz <- coxph(Surv(time, event) ~ rz, data = d); print(cox.zph(f_rz))

# ---- Nomogram (Figure 4) and calibration (Figure 5) --------------------------
nd <- d[complete.cases(d[, c("ridge_score", "stageS", "age", "hormone")]), c("time", "event", "ridge_score", "stageS", "age", "hormone")]
nd$stageS  <- factor(nd$stageS,  levels = c("Early (I-II)", "Advanced (III-IV)"))
nd$hormone <- factor(nd$hormone, levels = c("Non-functional", "Functional"))
label(nd$ridge_score) <- "Gene risk score"; label(nd$stageS) <- "ENSAT stage"
label(nd$age) <- "Age (years)"; label(nd$hormone) <- "Hormone status"
cat("Nomogram:", nrow(nd), "patients,", sum(nd$event), "events\n")       # expected 73 / 34
ddn <- datadist(nd); options(datadist = "ddn")
fn <- cph(Surv(time, event) ~ ridge_score + stageS + age + hormone, data = nd, x = TRUE, y = TRUE, surv = TRUE, time.inc = 365)
print(fn); print(cox.zph(coxph(Surv(time, event) ~ ridge_score + stageS + age + hormone, data = nd)))
set.seed(SEED); vn <- validate(fn, B = 500)
nomo_C <- c(apparent = vn["Dxy", "index.orig"] / 2 + .5, corrected = vn["Dxy", "index.corrected"] / 2 + .5)
print(round(nomo_C, 3))                                               # expected 0.799 / 0.779

sv <- Survival(fn)
nom <- nomogram(fn, age = seq(20, 80, by = 20), lp = FALSE,
                fun = list(function(x) sv(365, x), function(x) sv(1095, x)),
                funlabel = c("1-year RFS probability", "3-year RFS probability"), fun.at = c(.95, .9, .8, .7, .5, .3, .1))
png(file.path(DIR_FIG, "Figure4_nomogram.png"), width = 3000, height = 1800, res = 300)
plot(nom, xfrac = .28, cex.var = 1, cex.axis = .8, lmgp = .25); dev.off()

png(file.path(DIR_FIG, "Figure5_calibration.png"), width = 3200, height = 1600, res = 300)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2, 1))
for (u in c(365, 1095)) {
  fu <- cph(Surv(time, event) ~ ridge_score + stageS + age + hormone, data = nd, x = TRUE, y = TRUE, surv = TRUE, time.inc = u)
  set.seed(SEED); cal <- calibrate(fu, u = u, B = 200, cmethod = "KM", m = ifelse(u == 365, 18, 24))
  plot(cal, xlim = c(0, 1), ylim = c(0, 1), subtitles = FALSE,
       xlab = sprintf("Nomogram-predicted %d-year RFS", u / 365), ylab = sprintf("Observed %d-year RFS (Kaplan-Meier)", u / 365))
  abline(0, 1, lty = 2, col = "grey50")
  mtext(ifelse(u == 365, "A", "B"), side = 3, adj = -0.15, line = 0.5, font = 2, cex = 1.4)
  legend("bottomright", bty = "n", cex = .8, pch = c(19, 4), col = c("black", "blue"), legend = c("Apparent", "Bias-corrected (B=200)"))
}
dev.off()
if (requireNamespace("polspline", quietly = TRUE)) {
  for (u in c(365, 1095)) {
    fu <- cph(Surv(time, event) ~ ridge_score + stageS + age + hormone, data = nd, x = TRUE, y = TRUE, surv = TRUE, time.inc = u)
    set.seed(SEED); print(calibrate(fu, u = u, B = 200, cmethod = "hare"))   # mean |error| 0.063 / 0.076
  }
}

# ---- Table 2 -----------------------------------------------------------------
c_ci <- function(f) { cc <- concordance(f); se <- sqrt(cc$var); fmt_ci(cc$concordance, cc$concordance - 1.96 * se, cc$concordance + 1.96 * se) }
hr_ci <- function(f, i = 1) { s <- summary(f)$conf.int; fmt_ci(s[i, 1], s[i, 3], s[i, 4], 2) }
f_un <- coxph(S ~ predict(m6)); f_7 <- coxph(S ~ scale(score_z(d, b7)))
tab2 <- data.frame(
  Model = c("Six-gene ridge risk score (final model)", "Six-gene unpenalized Cox, apparent",
            "Six-gene unpenalized Cox, bootstrap-corrected", "Full pipeline, nested CV (mean; 2.5-97.5th pct)",
            "Seven-gene ridge (including CDK6)", "Clinical model (stage, age, hormone, resection)",
            "Clinical model + risk score", "Nomogram, apparent / bootstrap-corrected"),
  n_events = c("79 / 37", "79 / 37", "79 / 37", "79 / 37", "79 / 37",
               paste(nrow(dm), "/", sum(dm$event)), paste(nrow(dm), "/", sum(dm$event)), paste(nrow(nd), "/", sum(nd$event))),
  C = c(c_ci(f_rz), c_ci(f_un), sprintf("%.3f", boot_summary["C_corrected"]),
        fmt_ci(nested_summary["mean"], nested_summary["2.5%"], nested_summary["97.5%"]), c_ci(f_7),
        sprintf("%.3f", added["C_clin"]), sprintf("%.3f", added["C_full"]),
        sprintf("%.3f / %.3f", nomo_C["apparent"], nomo_C["corrected"])),
  HR_perSD = c(hr_ci(f_rz), "", "", "", hr_ci(f_7), "", hr_ci(full), ""),
  p = c(signif(summary(f_rz)$coefficients[5], 2), "", "", "", signif(summary(f_7)$coefficients[5], 2), "",
        signif(added["LRT_p"], 2), ""))
write.csv(tab2, file.path(DIR_RES, "Table2_TCGA_performance.csv"), row.names = FALSE)
print(tab2)
