###############################################################################
# 03_ppi_hub_selection.R
# STRING network of the 300 most significant DEGs, centrality metrics,
# hub-gene selection (Figure 2, Supplementary Table S2).
###############################################################################

sig_for_ppi <- sig_genes_annotated %>% distinct(SYMBOL, .keep_all = TRUE) %>% arrange(padj) %>% slice_head(n = 300)

string_db <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400, input_directory = DIR_STRING)  # files cached between runs
mapped <- string_db$map(sig_for_ppi, "SYMBOL", removeUnmappedRows = TRUE)
cat("Mapped to STRING:", nrow(mapped), "/ 300\n")                    # expected 261
ppi_network <- string_db$get_interactions(mapped$STRING_id)
id2sym <- setNames(mapped$SYMBOL, mapped$STRING_id)
ppi_network <- ppi_network %>% mutate(from_symbol = id2sym[from], to_symbol = id2sym[to]) %>%
  filter(!is.na(from_symbol), !is.na(to_symbol), from_symbol != to_symbol)
write.csv(ppi_network, file.path(DIR_RES, "PPI_network_edges.csv"), row.names = FALSE)

g <- simplify(graph_from_data_frame(ppi_network[, c("from_symbol", "to_symbol")], directed = FALSE))
cat("Network:", vcount(g), "nodes,", ecount(g), "unique edges\n")   # expected 183 / 365
centrality_df <- data.frame(gene = V(g)$name, degree = degree(g),
                            betweenness = betweenness(g, normalized = TRUE),
                            closeness = closeness(g, normalized = TRUE))

# ALB and mitochondrially encoded genes are frequent non-specific STRING hubs
centrality_filtered <- centrality_df %>% filter(!grepl("^ALB$|^MT-", gene)) %>% arrange(desc(degree), desc(betweenness))
centrality_filtered$rank <- seq_len(nrow(centrality_filtered))

# ---- Hub-gene selection ------------------------------------------------------
# Rule (formalized retrospectively): among the 20 highest-ranked nodes, retain genes
# annotated to GO:0051301 (cell division). TYMS (no cell-division annotation) was
# added on biological grounds (thymidylate synthesis for S-phase DNA replication).
top20 <- head(centrality_filtered$gene, 20)
cell_division <- unique(unlist(mget("GO:0051301", org.Hs.egGO2ALLEGS)))
top20_entrez <- AnnotationDbi::mapIds(org.Hs.eg.db, top20, "ENTREZID", "SYMBOL")
rule_genes <- top20[top20_entrez %in% cell_division]
selected_genes <- union(rule_genes, "TYMS")
cat("Rule-based genes:", paste(rule_genes, collapse = ", "), "\nSelected genes:", paste(selected_genes, collapse = ", "), "\n")
if (!setequal(selected_genes, GENES7)) warning("Selected genes differ from the published set: check STRING / GO versions.")

tabS2 <- head(centrality_filtered, 50) %>%
  left_join(sig_for_ppi[, c("SYMBOL", "log2FoldChange", "padj")], by = c("gene" = "SYMBOL")) %>%
  mutate(GO_0051301 = gene %in% AnnotationDbi::mapIds(org.Hs.eg.db, as.character(cell_division), "SYMBOL", "ENTREZID"),
         selected = gene %in% GENES7)
write.csv(tabS2, file.path(DIR_RES, "TableS2_PPI_centrality_top50.csv"), row.names = FALSE)

# ---- Figure 2: top-50 subnetwork ----------------------------------------------
g_sub <- induced_subgraph(g, vids = which(V(g)$name %in% head(centrality_filtered$gene, 50)))
set.seed(42)
p_net <- ggraph(g_sub, layout = "fr") +
  geom_edge_link(alpha = 0.3, colour = "grey60") +
  geom_node_point(aes(size = igraph::degree(g_sub)), colour = "#D55E00") +
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +
  theme_void() + labs(size = "Degree")
ggsave(file.path(DIR_FIG, "Figure2_PPI_top50.png"), p_net, width = 12, height = 10, dpi = 300)