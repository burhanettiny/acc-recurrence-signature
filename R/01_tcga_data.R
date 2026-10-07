###############################################################################
# 01_tcga_data.R
# TCGA-ACC RNA-seq download, clinical data (BCR Biotab), recurrence definition,
# recurrence-free survival (RFS) and clinical covariates.
###############################################################################

# ---- 1. RNA-seq (STAR-Counts) ------------------------------------------------
query_rna <- GDCquery(project = "TCGA-ACC",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")
GDCdownload(query_rna, directory = GDC_DIR)
acc_data   <- GDCprepare(query_rna, directory = GDC_DIR)
raw_counts <- assay(acc_data, "unstranded")
cat("RNA-seq:", nrow(raw_counts), "genes x", ncol(raw_counts), "samples\n")

# ---- 2. Clinical tables ------------------------------------------------------
# The BCR Biotab files used in the manuscript are bundled in data/raw_clinical_biotab
# (see MANIFEST.txt for GDC file identifiers and MD5 checksums).
patient_main <- read_biotab("patient_acc")
nte_main     <- read_biotab("nte_acc")
nte_v4       <- read_biotab("follow_up_v4.0_nte_acc")

# Follow-up (vital status, days to last follow-up / death) from the GDC clinical
# endpoint. A snapshot of the values used in the manuscript is bundled; it is used
# when present because GDC clinical records are updated over time.
snap_file <- file.path("data", "tcga_followup_snapshot.csv")
if (file.exists(snap_file)) {
  followup_time_df <- read.csv(snap_file, stringsAsFactors = FALSE)
} else {
  followup_time_df <- GDCquery_clinic(project = "TCGA-ACC", type = "clinical") %>%
    dplyr::select(submitter_id, vital_status, days_to_last_follow_up, days_to_death) %>%
    distinct()
}

# ---- 3. Recurrence groups (any recorded new tumor event, NTE) ----------------
all_patients <- unique(patient_main$bcr_patient_barcode[grepl("^TCGA-", patient_main$bcr_patient_barcode)])

extract_recurrent_patients <- function(df) {
  df <- df %>% filter(grepl("^TCGA-", bcr_patient_barcode))
  ind <- grep("new_tumor_event|new_neoplasm_event_type", colnames(df), value = TRUE, ignore.case = TRUE)
  if (length(ind) == 0) return(unique(df$bcr_patient_barcode))
  df %>% filter(if_any(all_of(ind), ~ toupper(trimws(.)) %in% "YES" |
                         (!is.na(.) & . != "" & toupper(.) != "NO"))) %>%
    pull(bcr_patient_barcode) %>% unique()
}
recurrent_patients <- unique(c(extract_recurrent_patients(nte_main), extract_recurrent_patients(nte_v4)))
clinical_grouped <- data.frame(submitter_id = all_patients, stringsAsFactors = FALSE) %>%
  mutate(recurrence_group = if_else(submitter_id %in% recurrent_patients, "Recurrent", "Non-recurrent"))
cat("Patients with clinical data:", length(all_patients), "| with NTE:", length(recurrent_patients), "\n")

# ---- 4. Match RNA-seq samples (primary tumours only) -------------------------
sample_info <- as.data.frame(colData(acc_data))
sample_info$submitter_id <- substr(sample_info$barcode, 1, 12)
sample_info <- sample_info %>%
  left_join(clinical_grouped, by = "submitter_id") %>%
  filter(!is.na(recurrence_group), shortLetterCode == "TP")
sample_info$recurrence_group <- factor(sample_info$recurrence_group, levels = c("Non-recurrent", "Recurrent"))
print(table(sample_info$recurrence_group))   # expected: 42 Non-recurrent, 37 Recurrent

# ---- 5. Recurrence-free survival ---------------------------------------------
extract_recurrence_time <- function(df) {
  df %>% filter(grepl("^TCGA-", bcr_patient_barcode)) %>%
    mutate(d = num(days_to_new_tumor_event_after_initial_treatment)) %>%
    filter(!is.na(d)) %>% group_by(bcr_patient_barcode) %>%
    summarise(days_to_recurrence = min(d), .groups = "drop")
}
recurrence_time_df <- bind_rows(extract_recurrence_time(nte_main), extract_recurrence_time(nte_v4)) %>%
  group_by(bcr_patient_barcode) %>% summarise(days_to_recurrence = min(days_to_recurrence), .groups = "drop")

survival_df <- data.frame(submitter_id = unique(sample_info$submitter_id), stringsAsFactors = FALSE) %>%
  left_join(recurrence_time_df, by = c("submitter_id" = "bcr_patient_barcode")) %>%
  left_join(followup_time_df, by = "submitter_id") %>%
  mutate(event = if_else(!is.na(days_to_recurrence), 1, 0),
         time = case_when(!is.na(days_to_recurrence) ~ days_to_recurrence,
                          !is.na(days_to_last_follow_up) ~ days_to_last_follow_up,
                          !is.na(days_to_death) ~ days_to_death, TRUE ~ NA_real_)) %>%
  filter(!is.na(time), time > 0)
cat("RFS dataset:", nrow(survival_df), "patients,", sum(survival_df$event), "events\n")

# ---- 6. Clinical covariates --------------------------------------------------
nte_first <- nte_main %>% filter(grepl("^TCGA-", bcr_patient_barcode)) %>%
  mutate(day = num(days_to_new_tumor_event_after_initial_treatment)) %>%
  arrange(bcr_patient_barcode, day) %>% distinct(bcr_patient_barcode, .keep_all = TRUE) %>%
  dplyr::select(bcr_patient_barcode, first_nte_type = new_tumor_event_type)

clin <- patient_main %>% filter(grepl("^TCGA-", bcr_patient_barcode)) %>%
  transmute(submitter_id = bcr_patient_barcode,
            age = num(age_at_initial_pathologic_diagnosis),
            sex = tolower(gender),
            stage_raw = ajcc_pathologic_tumor_stage,
            stageS = case_when(grepl("Stage I$|Stage II$", stage_raw) ~ "Early (I-II)",
                               grepl("Stage III|Stage IV", stage_raw) ~ "Advanced (III-IV)"),
            hormone = case_when(grepl("Not Available|Unknown", history_adrenal_hormone_excess, ignore.case = TRUE) |
                                  is.na(history_adrenal_hormone_excess) ~ NA_character_,
                                toupper(trimws(history_adrenal_hormone_excess)) == "NONE" ~ "Non-functional",
                                TRUE ~ "Functional"),
            residual_tumor = ifelse(residual_tumor %in% c("R0", "R1", "R2", "RX"), residual_tumor, NA),
            nonR0 = case_when(residual_tumor == "R0" ~ 0, residual_tumor %in% c("R1", "R2") ~ 1),
            mitotane_adjuvant = ifelse(pharm_tx_mitotane_adjuvant %in% c("YES", "NO"), pharm_tx_mitotane_adjuvant, NA)) %>%
  left_join(nte_first, by = c("submitter_id" = "bcr_patient_barcode"))

write.csv(survival_df, file.path(DIR_RES, "TCGA_survival_dataset.csv"), row.names = FALSE)
