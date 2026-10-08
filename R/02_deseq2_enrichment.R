###############################################################################
# 02_deseq2_enrichment.R
# Differential expression (recurrent vs non-recurrent), stage-adjusted DESeq2,
# GO / KEGG enrichment (Figure 1, Supplementary Table S1).
###############################################################################

counts_matched <- raw_counts[, sample_info$barcode]
stopifnot(all(colnames(counts_matched) == sample_info$barcode))

dds <- DESeqDataSetFromMatrix(countData = round(counts_matched), colData = sample_info,
                              design = ~ recurrence_group)
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]
cat("Genes after filtering:", nrow(dds), "\n")                       # expected 22,737
dds <- DESeq(dds)
res_df <- as.data.frame(results(dds, contrast = c("recurrence_group", "Recurrent", "Non-recurrent"))) %>%
  rownames_to_column("gene_id") %>% arrange(padj)
res_df$gene_symbol <- rowData(acc_data)$gene_name[match(res_df$gene_id, rownames(acc_data))]
write.csv(res_df, file.path(DIR_RES, "DEG_recurrence_vs_nonrecurrence.csv"), row.names = FALSE)

# Expression matrices: VST (blind = FALSE) for survival modelling of the hub genes;
# blind VST for the genome-wide nested cross-validation and the BUB1B-PINK1 score.
vst_matrix   <- assay(vst(dds, blind = FALSE))
vst_blind    <- assay(vst(dds, blind = TRUE))
symbol_of    <- setNames(rowData(acc_data)$gene_name, rownames(acc_data))

# ---- Significant DEGs --------------------------------------------------------
sig_genes <- res_df %>% filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1)
cat("Significant DEGs:", nrow(sig_genes), "| up:", sum(sig_genes$log2FoldChange > 0),
    "| down:", sum(sig_genes$log2FoldChange < 0), "\n")           # expected 2,101 (1,072 / 1,029)
sig_genes$ensembl_clean <- sub("\\..*$", "", sig_genes$gene_id)
gene_conversion <- bitr(sig_genes$ensembl_clean, fromType = "ENSEMBL",
                        toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
sig_genes_annotated <- sig_genes %>% dplyr::select(-gene_symbol) %>%
  left_join(gene_conversion, by = c("ensembl_clean" = "ENSEMBL")) %>% filter(!is.na(ENTREZID))
write.csv(sig_genes_annotated, file.path(DIR_RES, "DEG_significant_annotated.csv"), row.names = FALSE)

# ---- Figure 1: volcano plot --------------------------------------------------
if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  rv  <- res_df[!is.na(res_df$padj), ]
  sg  <- rv[rv$padj < 0.05 & abs(rv$log2FoldChange) > 1, ]
  lab <- unique(c(head(sg$gene_symbol[order(sg$padj)], 15), GENES6))   # top 15 DEGs + signature genes
  p_vol <- EnhancedVolcano::EnhancedVolcano(rv, lab = rv$gene_symbol, x = "log2FoldChange", y = "padj",
             selectLab = lab, pCutoff = 0.05, FCcutoff = 1, labSize = 3.5,
             drawConnectors = TRUE, widthConnectors = 0.4, max.overlaps = Inf, boxedLabels = FALSE,
             ylab = bquote(~-Log[10] ~ "adjusted" ~ italic(P)),
             legendLabels = c("NS", expression(Log[2] ~ FC), "adj. P", expression(adj. ~ P ~ and ~ log[2] ~ FC)),
             title = "Recurrent vs non-recurrent ACC", subtitle = "TCGA-ACC cohort, n = 79", caption = "")
  ggsave(file.path(DIR_FIG, "Figure1_volcano.png"), p_vol, width = 9, height = 7.5, dpi = 300)
}

# ---- GO / KEGG (Supplementary Table S1) -------------------------------------
entrez <- unique(sig_genes_annotated$ENTREZID)
go_all <- enrichGO(entrez, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05, readable = TRUE)
kegg <- setReadable(enrichKEGG(entrez, organism = "hsa", pAdjustMethod = "BH",
                               pvalueCutoff = 0.05, qvalueCutoff = 0.05), OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
cat("GO terms:", nrow(as.data.frame(go_all)), "| KEGG pathways:", nrow(as.data.frame(kegg)), "\n")  # expected 664 / 39 (1,829 DEGs mapped to Entrez IDs)
write.csv(as.data.frame(go_all), file.path(DIR_RES, "TableS1a_GO_full.csv"), row.names = FALSE)
write.csv(as.data.frame(kegg),   file.path(DIR_RES, "TableS1b_KEGG_full.csv"), row.names = FALSE)

# ---- Stage-adjusted DESeq2 (differential expression only) -------------------
si_st <- sample_info %>% left_join(clin[, c("submitter_id", "stageS")], by = "submitter_id") %>% filter(!is.na(stageS))
si_st$stage_simplified <- factor(si_st$stageS, levels = c("Early (I-II)", "Advanced (III-IV)"))
dds_st <- DESeqDataSetFromMatrix(round(raw_counts[, si_st$barcode]), si_st, ~ stage_simplified + recurrence_group)
dds_st <- DESeq(dds_st[rowSums(counts(dds_st) >= 10) >= 10, ])
res_st <- as.data.frame(results(dds_st, contrast = c("recurrence_group", "Recurrent", "Non-recurrent"))) %>%
  rownames_to_column("gene_id")
res_st$gene_symbol <- symbol_of[res_st$gene_id]
print(res_st[res_st$gene_symbol %in% GENES6, c("gene_symbol", "log2FoldChange", "padj")])
write.csv(res_st, file.path(DIR_RES, "DEG_recurrence_stage_adjusted.csv"), row.names = FALSE)
