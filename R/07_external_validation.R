###############################################################################
# 07_external_validation.R
# GEO cohorts GSE10927, GSE19750, GSE76021: preprocessing, primary-carcinoma
# definitions, Table 3, Supplementary Tables S4/S5/S8/S9, head-to-head comparison
# with BUB1B-PINK1, Figure 6.
###############################################################################

genes_all <- c(GENES7, "BUB1B", "PINK1")
geo <- function(acc) getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE)[[1]]

# ---- GSE10927 (HG-U133 Plus 2.0; log10 MAS5; OS) -----------------------------
g10 <- geo("GSE10927"); p10 <- pData(g10)
s10 <- p10 %>% filter(`tumor stage:ch1` != "not applicable", !is.na(`tumor stage:ch1`)) %>%
  transmute(gsm = geo_accession, time = num(`years to last followup:ch1`) * 365.25,
            event = case_when(toupper(trimws(`dead or alive at last followup:ch1`)) == "DEAD" ~ 1,
                              toupper(trimws(`dead or alive at last followup:ch1`)) == "ALIVE" ~ 0),
            side = `side of body:ch1`, clin_note = `cinical characteristics:ch1`) %>%
  filter(!is.na(time), !is.na(event), time > 0)
e10 <- gene_matrix(exprs(g10), hgu133plus2.db, genes_all, s10$gsm)
D10 <- cbind(s10, e10)
D10$probably_nonprimary <- D10$side == "not applicable" |
  grepl("History of Adrenocortical carcinoma|Status post adrenalectomy|adrenalectomy 4 years prior|Liver tumor|hepatic lobe",
        D10$clin_note, ignore.case = TRUE)

# ---- GSE19750 (HG-U133 Plus 2.0; linear MAS5 -> log2; OS) --------------------
g19 <- geo("GSE19750"); p19 <- pData(g19)
s19 <- p19 %>% filter(`Stage:ch1` != "NA", !is.na(`Stage:ch1`)) %>%
  transmute(gsm = geo_accession, stage = `Stage:ch1`, time = num(`survival in years:ch1`) * 365.25,
            event = case_when(toupper(trimws(`survival status:ch1`)) == "DEAD" ~ 1,
                              toupper(trimws(`survival status:ch1`)) == "ALIVE" ~ 0)) %>%
  filter(!is.na(time), !is.na(event), time > 0)
e19 <- log2(gene_matrix(exprs(g19), hgu133plus2.db, genes_all, s19$gsm))
D19 <- cbind(s19, e19); D19$primary <- !grepl("Recurrence|Metasta", D19$stage)

# ---- GSE76021 (HG-U133A; log2; pediatric; EFS) -------------------------------
g76 <- geo("GSE76021"); p76 <- pData(g76)
s76 <- data.frame(gsm = p76$geo_accession, histology = p76$`histology:ch1`,
                  stage = sub("Stage: ", "", p76$characteristics_ch1.1),
                  time = num(p76$`efs.time:ch1`), event = num(p76$`efs.event:ch1`))
e76 <- gene_matrix(exprs(g76), hgu133a.db, genes_all, s76$gsm)
D76 <- cbind(s76, e76)

cat("GSE10927:", nrow(D10), "| GSE19750:", nrow(D19), "(primary", sum(D19$primary), ") | GSE76021:",
    nrow(D76), "(ACC", sum(D76$histology == "ACC"), ")\n")    # expected 24 | 21 (13) | 29 (19)

tabS8 <- data.frame(Cohort = c("GSE10927", "GSE19750", "GSE76021"),
                    Platform = c("GPL570", "GPL570", "GPL96"),
                    Arrays = c(ncol(g10), ncol(g19), ncol(g76)),
                    With_survival = c(nrow(D10), nrow(D19), nrow(D76)),
                    Primary_analysis_set = c(nrow(D10), sum(D19$primary), sum(D76$histology == "ACC")),
                    Endpoint = c("OS", "OS", "EFS"))
write.csv(tabS8, file.path(DIR_RES, "TableS8_external_cohorts.csv"), row.names = FALSE)

# ---- Primary analysis sets and sensitivity sets -----------------------------
TCGA <- list(D = d, t = d$time, e = d$event)
mk <- function(D) list(D = D, t = D$time, e = D$event)
coh_main <- list(TCGA = TCGA, GSE10927 = mk(D10), GSE19750 = mk(D19[D19$primary, ]), GSE76021 = mk(D76[D76$histology == "ACC", ]))
coh_sens <- list("GSE10927 probable primary" = mk(D10[!D10$probably_nonprimary, ]),
                 "GSE19750 all samples" = mk(D19), "GSE76021 all tumours" = mk(D76))

w6o <- coefs_6gene * apply(as.matrix(d[, GENES6]), 2, sd)   # original weights, per SD
models <- list(list(b6, "Final six-gene ridge"), list(w6o, "Six-gene, original weights"),
               list(setNames(rep(1, 6), GENES6), "Six-gene, unweighted"),
               list(b7, "Seven-gene ridge (incl. CDK6)"), list(b5, "GO:0051301 rule, CDK6 excluded"))

tab3 <- rbind(do.call(rbind, lapply(names(coh_main), function(nm) { co <- coh_main[[nm]]
                 cox_row(co$t, co$e, score_z(co$D, b6), "Final six-gene ridge", nm) })),
              do.call(rbind, lapply(names(coh_sens), function(nm) { co <- coh_sens[[nm]]
                 cox_row(co$t, co$e, score_z(co$D, b6), "Final six-gene ridge (sensitivity)", nm) })))
print(tab3); write.csv(tab3, file.path(DIR_RES, "Table3_external.csv"), row.names = FALSE)

tabS459 <- do.call(rbind, lapply(models, function(m) do.call(rbind, lapply(names(coh_main), function(nm) {
  co <- coh_main[[nm]]; cox_row(co$t, co$e, score_z(co$D, m[[1]]), m[[2]], nm) }))))
write.csv(tabS459, file.path(DIR_RES, "TableS4_S5_S9_models.csv"), row.names = FALSE)

# Stage-adjusted model in the pediatric carcinomas
acc76 <- coh_main$GSE76021$D; acc76$rz <- as.vector(scale(score_z(acc76, b6)))
acc76$adv <- as.numeric(acc76$stage %in% c("III", "IV"))
print(summary(coxph(Surv(time, event) ~ rz + adv, data = acc76))$conf.int)

# ---- Head-to-head comparison with BUB1B-PINK1 --------------------------------
h2h <- do.call(rbind, lapply(names(coh_main)[-1], function(nm) {
  co <- coh_main[[nm]]; y <- Surv(co$t, co$e)
  s <- score_z(co$D, b6); bp <- co$D$BUB1B - co$D$PINK1
  fs <- coxph(y ~ s); fb <- coxph(y ~ bp); fsb <- coxph(y ~ s + bp)
  data.frame(Cohort = nm, C_score = concordance(fs)$concordance, C_BUB1B_PINK1 = concordance(fb)$concordance,
             p_score_added = anova(fb, fsb)[2, 4], p_BUB1B_PINK1_added = anova(fs, fsb)[2, 4])
}))
print(h2h); write.csv(h2h, file.path(DIR_RES, "TableS6c_BUB1B_PINK1_external.csv"), row.names = FALSE)

# ---- Figure 6 ------------------------------------------------------------------
km_panel <- function(t_years, e, score, title, endpoint, pc = c(0, 0.06)) {
  df <- data.frame(t = t_years, e = e, grp = factor(ifelse(score > median(score), "High risk", "Low risk"),
                                                     levels = c("Low risk", "High risk")))
  ggsurvplot(survfit(Surv(t, e) ~ grp, data = df), data = df, pval = TRUE, pval.size = 3.5, pval.coord = pc,
             risk.table = TRUE, risk.table.height = 0.25, risk.table.title = "Number at risk",
             tables.y.text = FALSE, tables.col = "strata", fontsize = 3,
             tables.theme = theme_cleantable() + theme(plot.title = element_text(size = 9)),
             palette = c("#2E6DB4", "#C0392B"), legend.labs = c("Low risk", "High risk"), legend.title = "",
             xlab = "Time (years)", ylab = paste(endpoint, "probability"), title = title, break.time.by = 2,
             censor.size = 3, font.title = c(11, "bold"), font.x = 10, font.y = 10, font.tickslab = 9, font.legend = 9)
}
m <- coh_main
p6 <- list(km_panel(m$TCGA$t / 365.25, m$TCGA$e, score_z(d, b6), "A  TCGA-ACC (discovery, n=79)", "RFS"),
           km_panel(m$GSE10927$t / 365.25, m$GSE10927$e, score_z(m$GSE10927$D, b6), "B  GSE10927 (n=24)", "OS"),
           km_panel(m$GSE19750$t / 365.25, m$GSE19750$e, score_z(m$GSE19750$D, b6), "C  GSE19750, primary ACC (n=13)", "OS"),
           km_panel(m$GSE76021$t, m$GSE76021$e, score_z(m$GSE76021$D, b6), "D  GSE76021, pediatric ACC (n=19)", "EFS", pc = c(10, 0.95)))
png(file.path(DIR_FIG, "Figure6_KM_external.png"), width = 3600, height = 3200, res = 300)
print(arrange_ggsurvplots(p6[c(1, 3, 2, 4)], ncol = 2, nrow = 2, print = FALSE)); dev.off()
