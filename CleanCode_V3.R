#Load libraries
library(Seurat)
library(dplyr)
library(presto)
library(ggplot2)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(stringr)
library(harmony)
library(msigdbr)
library(fgsea)
library(homologene)
library(UpSetR)
library(VennDiagram)
library(ComplexUpset)
library(enrichR)
library(tibble)
library(GO.db)
library(AnnotationDbi)
library(scales)
library(tidyr)
library(dorothea)
library(viper)
library(purrr)
library(pheatmap)

# Figure 1 Analysis
###############################################################################
AT2_1 <- readRDS("LTRC_AT2only.rds")
AT2_2 <- readRDS("Penn_AT2only.rds")

# Figure 1A: LTRC Volcano plot
Idents(AT2_1) <- "disease_state"

DEG <- FindMarkers(
  AT2_1,
  ident.1 = "AATD",
  ident.2 = "COPD",
  only.pos = FALSE,
  min.pct = 0.25,
  test.use = "wilcox", 
  recorrect_umi = FALSE
)
write.csv(DEG, "LTRC_AATDvsCOPD_DEG.csv", row.names = TRUE)

sum(DEG$p_val_adj < 0.05 & DEG$avg_log2FC < -0.25, na.rm = TRUE) # < -.25 for COPD

genes_lable <- c("FOSB", "CXCL2", "HOPX", "FAM13A", "ZFAND5", "XBP1", "NFKBIZ", "ICAM1", "NFKB1", "SOD2", "MAPKAPK2", "CFLAR", "KLF6", "JUND", "CTNNB1", "FOS", "EIF4A1", "SFTPA1")

volcano <- EnhancedVolcano(
  DEG,
  lab = rownames(DEG),
  x = "avg_log2FC",
  y = "p_val_adj",
  selectLab = genes_lable, 
  pCutoff = 0.05,
  FCcutoff = 0.25,
  pointSize = 1.5,
  labSize = 5.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.8,
  legendPosition = 'right',
  legendLabSize = 8,
  legendIconSize = 3.0,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  ylab = bquote(~-Log[10]~ 'padj'),
  xlab = bquote(~Log[2]~ 'FC'),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey30',
  max.overlaps = 105, 
) + 
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
volcano

ggsave(
  filename = "LTRC_AT2_Volcano_AATDvsCOPD.pdf",
  plot = volcano,
  width = 10, height = 6)

#Figure 1B: LTRC Hallmark GSEA 
msig.H <- msigdbr(species = "Homo sapiens", category = "H")
msig.H <- split(msig.H$gene_symbol, msig.H$gs_name)

gene_ranks <- DEG$avg_log2FC
names(gene_ranks) <- rownames(DEG)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

fgsea_res <- fgsea(
  pathways = msig.H,
  stats = gene_ranks,
  minSize = 5,
  maxSize = 500,    scoreType = "std"
)
fgsea_res_df <- fgsea_res %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))

write.csv(fgsea_res_df, "LTRC_Hallmark.csv", row.names = FALSE)

fgsea_res$Significance <- factor(
  ifelse(fgsea_res$padj < 0.05, "padj < 0.05", "Not sig."),
  levels = c("padj < 0.05", "Not sig.")
)
fgsea_res$`-log(P-value)` <- -log10(fgsea_res$pval)
fgsea_res$pretty_pathway <- gsub("HALLMARK_", "", fgsea_res$pathway)
fgsea_res$pretty_pathway <- gsub("_", " ", fgsea_res$pretty_pathway)
fgsea_res_ord <- fgsea_res[order(abs(fgsea_res$NES), decreasing = TRUE), ]
fgseaRes_top <- head(fgsea_res_ord, 20)

p <- ggplot(
  fgseaRes_top,
  aes(
    x = reorder(str_wrap(str_sub(pretty_pathway, 1, 80), 40), NES),
    y = NES,
    fill = Significance
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  ggtitle("HALLMARK GSEA (AATD vs COPD)") +
  theme(
    axis.text = element_text(size = 7, face = "bold"),
    axis.title.y = element_blank()
  ) +
  scale_fill_manual(
    name = "Significance",
    values = c("padj < 0.05" = "blue", "Not sig." = "grey")
  ) +
  ylab("Normalized Enrichment Score (NES)")
p

ggsave(
  filename = "LTRC_AT2_Hallmark_AATDvsCOPD.pdf",  
  plot = p,
  width = 6, height = 7
)

#Figure 1C: Regulon Analysis
dorothea_regulon_human <- get(data("dorothea_hs", package = "dorothea"))

regulon <- dorothea_regulon_human %>%
  dplyr::filter(tf %in% c("AR","ARNTL","ATF1","ATF2","ATF4","ATF6","CDX2","CEBPA","CEBPB","CREB1","CTCF","E2F1","E2F2","E2F3","E2F4","EGR1","ELK1","EPAS1","ERG","ESR1","ESR2","ETS1","ETS2","ETV4","FLI1","FOS","FOSL1","FOSL2","FOXA1","FOXL2","FOXM1","FOXO1","FOXO3","FOXO4","GATA1","GATA2","GATA3","GLI2","HIF1A","HNF1A","HNF4A","IRF1","JUN","JUND","KLF4","KMT2A","LEF1","MITF","MYB","MYC","MYCN","NFATC2","NFIC","NFKB1","NR2F2","NR3C1","NR5A1","PAX6","PAX8","PGR","POU2F1","PPARA","PPARG","RARA","REL","RELA","RFX5","RUNX1","RXRA","SMAD3","SMAD4","SOX10","SOX2","SOX9","SP1","SP3","SPI1","SREBF1","SREBF2","SRF","STAT1","STAT3","STAT5A","STAT5B","STAT6","TAL1","TCF7L2","TFAP2A","TFAP2C","TP53","TWIST1","USF1","USF2","VDR","WT1","YY1"))

regulon_viper <- regulon %>%
  split(.$tf) %>%
  map(function(df) {
    list(
      tfmode = setNames(df$mor, df$target),
      likelihood = setNames(rep(1, nrow(df)), df$target)
    )
  })

Idents(AT2_1) <- "disease_state"
AT2_1 <- ScaleData(AT2_1)

expr_mat <- GetAssayData(
  AT2_1,
  assay = "SCT",
  layer = "scale.data"
)

viper_res <- viper(
  expr_mat,
  regulon_viper,
  method = "scale",
  minsize = 5,
  verbose = FALSE
)

AT2_1[["VIPER"]] <- CreateAssayObject(data = viper_res)
DefaultAssay(AT2_1) <- "VIPER"

AT2_1 <- ScaleData(AT2_1)

viper_scores_df <- GetAssayData(
  AT2_1,
  assay = "VIPER",
  layer = "scale.data"     
) %>%
  as.matrix() %>%     
  t()                
head(viper_scores_df)

viper_scores_df_fixed <- viper_scores_df
rownames(viper_scores_df_fixed) <- sub("-", ".", rownames(viper_scores_df_fixed))

viper_scores_clusters <- viper_scores_df_fixed %>%
  data.frame(check.names = FALSE) %>%
  rownames_to_column("cell") %>%
  pivot_longer(cols = -cell, names_to = "tf", values_to = "activity") %>%
  inner_join(CellsClusters, by = "cell")
head(viper_scores_clusters)

CellsClusters <- data.frame(cell = names(Idents(AT2_1)),
                            cell_type = as.character(AT2_1$disease_state),
                            stringsAsFactors = FALSE,
                            check.names=FALSE) %>% dplyr::mutate(., cell = sapply(.$cell, function(x) sub("-", ".", x)) %>% unname())
head(CellsClusters)

summarized_viper_scores <- viper_scores_clusters %>% 
  group_by(tf, cell_type) %>%
  summarise(avg = mean(activity),
            std = sd(activity))
head(summarized_viper_scores)

highly_variable_tfs <- summarized_viper_scores %>%
  group_by(tf) %>%
  mutate(var = var(avg))  %>%
  ungroup() %>%
  top_n(200, var) %>%
  distinct(tf)
head(highly_variable_tfs)

summarized_viper_scores_df <- summarized_viper_scores %>%
  semi_join(highly_variable_tfs, by = "tf") %>%
  dplyr::select(-std) %>%   
  spread(tf, avg) %>%
  data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
head(summarized_viper_scores_df)

mean_diff <- summarized_viper_scores_df["AATD", ] - summarized_viper_scores_df["COPD", ]

mean_diff_vec <- as.numeric(mean_diff)
names(mean_diff_vec) <- colnames(summarized_viper_scores_df)

up_tfs <- names(mean_diff_vec[mean_diff_vec > 0])
down_tfs <- names(mean_diff_vec[mean_diff_vec < 0])

ordered_tfs <- c(up_tfs, down_tfs)
viper_matrix_ordered <- summarized_viper_scores_df[, ordered_tfs]

palette_length <- 100
my_color <- colorRampPalette(c("#0517BA", "white", "#BA1505"))(palette_length)
my_breaks <- c(seq(min(viper_matrix_ordered), 0, length.out = ceiling(palette_length/2) + 1),
               seq(max(viper_matrix_ordered)/palette_length, max(viper_matrix_ordered), length.out = floor(palette_length/2)))

viper_hmap <- pheatmap(t(viper_matrix_ordered),
         fontsize=14,
         fontsize_row=8,
         color=my_color,
         breaks=my_breaks,
         main="AT2 Regulon Activity",
         angle_col=0,
         border_color=NA,
         cluster_rows=TRUE,  
         cluster_cols=TRUE)   
viper_hmap

ggsave(
  filename = "LTRC_AT2_Regulon_AATDvsCOPD.pdf",  
  plot = viper_hmap,
  width = 6, height = 8.5
)

# Supplement 1A: Penn Volcano plot 
Idents(AT2_2) <- "annot2"

DEG <- FindMarkers(
  AT2_2,
  ident.1 = "AATD",
  ident.2 = "COPD",
  only.pos = FALSE,
  min.pct = 0.25,
  test.use = "wilcox", 
  recorrect_umi = FALSE
)
write.csv(DEG, "Penn_AATDvsCOPD_DEG.csv", row.names = TRUE)

sum(DEG$p_val_adj < 0.05 & DEG$avg_log2FC < -0.25,na.rm = TRUE) # < -.25 for COPD

genes_lable <- c("FOSB", "CXCL2", "HOPX", "FAM13A", "ZFAND5", "XBP1", "NFKBIZ", "ICAM1", "NFKB1", "SOD2", "MAPKAPK2", "CFLAR", "KLF6", "JUND", "CTNNB1", "FOS", "EIF4A1", "SFTPA1")

volcano <- EnhancedVolcano(
  DEG,
  lab = rownames(DEG),
  x = "avg_log2FC",
  y = "p_val_adj",
  selectLab = genes_lable, 
  pCutoff = 0.05,
  FCcutoff = 0.25,
  pointSize = 1.5,
  labSize = 5.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.8,
  legendPosition = 'right',
  legendLabSize = 8,
  legendIconSize = 3.0,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  ylab = bquote(~-Log[10]~ 'padj'),
  xlab = bquote(~Log[2]~ 'FC'),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey30',
  max.overlaps = 105, 
) + 
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
volcano

ggsave(
  filename = "Penn_AT2_Volcano_AATDvsCOPD.pdf",
  plot = volcano,
  width = 10, height = 6)


#Figure 1D: Venn Diagram & EnrichR Penn vs LTRC
deg_ltrc <- read.csv("LTRC_AATDvsCOPD_DEG.csv", row.names = 1) 
deg_penn <- read.csv("Penn_AATDvsCOPD_DEG.csv", row.names = 1) 

genes_ltrc <- rownames( deg_ltrc[deg_ltrc$avg_log2FC < -0.25 & deg_ltrc$p_val_adj < 0.05, ] ) 
genes_penn <- rownames( deg_penn[deg_penn$avg_log2FC < -0.25 & deg_penn$p_val_adj < 0.05, ] ) 
gene_lists <- list( LTRC = genes_ltrc, Penn = genes_penn ) 
pdf("Penn_vs_LTRC_Venn.pdf", width = 6, height = 6) 
venn.plot <- venn.diagram( x = list( LTRC = genes_ltrc, Penn = genes_penn ), 
                           filename = NULL, fill = c("lightblue", "lightgreen"), 
                           alpha = 0.7, cex = 1.5, cat.cex = 0, main = "Shared Upregulated Genes\n(avg_log2FC > 0.25, padj < 0.05)" ) 
grid.newpage() 
grid.draw(venn.plot) 
dev.off()

shared_genes <- intersect(genes_ltrc, genes_penn) 
length(shared_genes) 
shared_genes

enrichr_results <- enrichr(
  genes = shared_genes,
  databases = "MSigDB_Hallmark_2020")
hallmark_df <- enrichr_results[["MSigDB_Hallmark_2020"]]
sig <- hallmark_df %>%
  filter(Adjusted.P.value < 0.05)
head(sig[, c("Term", "Adjusted.P.value", "Odds.Ratio", "Genes")], 10)
##############################################################################
#Figure 3 Analysis
##################################################################################
AT2 <- readRDS("iAT2_nointestinal.rds")

# Figure 3A: UMAP iAT2
umap <- DimPlot(
  object = AT2,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 0.5  
) +
  scale_color_manual(
    values = c(
      "PiZZ_Air" = "#F525E4",
      "PiMM_Air" = "blue"
    )
  ) +
  theme_classic() +     
  theme(
    plot.title = element_blank(),           
    axis.text = element_blank(),            
    axis.ticks = element_blank(),           
    legend.title = element_blank()          
  )
umap

ggsave(
  filename = "iAT2_UMAP.pdf",  
  plot = umap,
  width = 6, height = 5
)

#Figure 3B: Louvain Clusters
umap <- DimPlot(
  object = AT2,
  reduction = "umap",
  group.by = "SCT_snn_res.0.2",
  pt.size = 0.5  
) +
  theme_classic() +     
  theme(
    plot.title = element_blank(),           
    axis.text = element_blank(),            
    axis.ticks = element_blank(),           
    legend.title = element_blank()          
  )
umap

ggsave(
  filename = "iAT2_UMAP_clusters.pdf",  
  plot = umap,
  width = 6, height = 5)

#Figure 3C: AT2 & non-AT2 markers 
Idents(AT2) <- "SCT_snn_res.0.2"
genes <- c("NKX2-1","SFTPB","SLC34A2","ABCA3","CDX2","AFP", "TP63", "SCGB3A2","FOXJ1")

p <- DotPlot(
  object = AT2,
  features = genes
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme_classic() +
  labs(
    x = "SCT_snn_res.0.2 cluster",
    y = "Gene"
  )
p

ggsave(
  filename = "iAT2_markers_clusters.pdf",
  plot = p,
  width = 8,
  height = 4
)

#Figure 3D: Heatmap
Idents(AT2) <- "SCT_snn_res.0.2"
DEG <- FindAllMarkers(
  AT2,
  only.pos = TRUE,
  min.pct = 0.25,
  test.use = "wilcox"
)
write.csv(DEG, "iAT2_Clusters.csv")

top50_df <- DEG %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 50,
    with_ties = FALSE
  ) %>%
  ungroup()
top50_genes <- unique(top50_df$gene)

AT2 <- ScaleData(
  AT2,
  features = top50_genes,
  verbose = TRUE
)

g <- DoHeatmap(
  AT2,
  features = top50_genes
) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(
    axis.text.y = element_text(size = 4)
  )
g

ggsave(
  filename = "iAT2_heatmap_clusters.pdf",
  plot = g,
  width = 5,
  height = 8
)

#Figure 3E: Cluster Defining Genes UMAP
umap <- FeaturePlot(AT2, feature = "TOP2A")
umap <- FeaturePlot(AT2, feature = "HSPA6")
umap <- FeaturePlot(AT2, feature = "IFI44L")
umap <- FeaturePlot(AT2, feature = "TSPO")

ggsave(
  filename = "TSPO_UMAP.pdf",  
  plot = umap,
  width = 6, height = 5
)

#Figure 3F: Hallmark Pathways Defining Each Cluster
hallmark_df <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_sets <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

fgsea_results <- DEG %>%
  group_by(cluster) %>%
  group_modify(~ {
    
    # ranked gene vector for this cluster
    ranks <- .x %>%
      distinct(gene, avg_log2FC) %>%
      arrange(desc(avg_log2FC)) %>%
      deframe()
    
    fgsea(
      pathways = hallmark_sets,
      stats = ranks,
      minSize = 15,
      maxSize = 500,
      nperm = 10000
    )
  }) %>%
  ungroup()

top_hallmarks <- fgsea_results %>%
  filter(padj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = NES,
    n = 5,
    with_ties = FALSE
  ) %>%
  ungroup()
top_hallmarks %>%
  dplyr::select(cluster, pathway, NES, padj)
print(top_hallmarks, n=50)

#Figure 3G: Heatmap directly comparing clusters 4 & 1
AT2_subset <- subset(AT2, idents = c("4", "1"))

DEG <- FindMarkers(
  AT2,
  ident.1 = "4",
  ident.2 = "1",
  only.pos = FALSE,
  min.pct = 0.25,
  test.use = "wilcox"
)
write.csv(DEG, "iAT2_4vs1.csv")

top20_clust4 <- DEG %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::slice_head(n = 20) %>%
  rownames()

top20_clust1 <- DEG %>%
  dplyr::filter(avg_log2FC < 0) %>%
  dplyr::arrange(avg_log2FC) %>%   # most negative first
  dplyr::slice_head(n = 20) %>%
  rownames()

heatmap_genes <- unique(c(top20_clust4, top20_clust1))
length(heatmap_genes)  # should be 40

AT2_subset <- ScaleData(
  AT2_subset,
  features = heatmap_genes,
  verbose = FALSE
)

g <- DoHeatmap(
  AT2_subset,
  features = heatmap_genes,
) +
  scale_fill_gradientn(colors = c("blue", "white", "red")) +
  theme(
    axis.text.y = element_text(size = 6)
  )
g

ggsave(
  filename = "Cluster4vs1_heatmap.pdf",  
  plot = g,
  width = 5, height = 10
)

#Figure 3H: ER stress specific GO pathways cluster 4 vs 1
go <- msigdbr(
  species = "Homo sapiens",
  category = "C5",
) %>%
  dplyr::select(gs_name, gene_symbol)

go_terms_of_interest <- c("GOMF_PROTEIN_FOLDING_CHAPERONE", "GOBP_CHAPERONE_MEDIATED_AUTOPHAGY", "GOBP_REGULATION_OF_PROTEIN_FOLDING", "GOBP_CELLULAR_RESPONSE_TO_UNFOLDED_PROTEIN", "GOBP_RESPONSE_TO_ENDOPLASMIC_RETICULUM_STRESS", "GOBP_PROTEIN_FOLDING_IN_ENDOPLASMIC_RETICULUM", "GOBP_POSITIVE_REGULATION_OF_AUTOPHAGY", "GOBP_ERAD_PATHWAY",
                          "GOBP_REGULATION_OF_ERAD_PATHWAY", "GOBP_ATF6_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", "GOBP_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", "GOBP_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", "GOBP_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", "GOBP_REGULATION_OF_ENDOPLASMIC_RETICULUM_UNFOLDED_PROTEIN_RESPONSE", "GOBP_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", "GOBP_POSITIVE_REGULATION_OF_ERAD_PATHWAY", 
                          "GOBP_REGULATION_OF_PERK_MEDIATED_UNFOLDED_PROTEIN_RESPONSE", "GOBP_POSITIVE_REGULATION_OF_IRE1_MEDIATED_UNFOLDED_PROTEIN_RESPONSE")

go_sets <- go %>%
  dplyr::filter(gs_name %in% go_terms_of_interest)

DEG <- DEG %>% 
  tibble::rownames_to_column(var = "gene")
ranks_vec <- DEG %>%
  distinct(gene, avg_log2FC) %>%  
  arrange(desc(avg_log2FC))       
ranks_vec <- setNames(ranks_vec$avg_log2FC, ranks_vec$gene)

go_pathways <- split(go_sets$gene_symbol, go_sets$gs_name)

fgsea_results <- fgseaMultilevel(
  pathways = go_pathways,
  stats = ranks_vec,
  minSize = 10,
  maxSize = 500
)

fgsea_sig <- fgsea_results %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES))

print(fgsea_sig[, c("pathway", "padj", "NES", "leadingEdge")], n=25)

#Figure 3I: PERK Pathway Module Score Clusters + Supplement Other UPR Arms 
PERK <- c("ASNS", "ATF3", "ATF4", "ATF6", "CARS", "CEBPB", "DDIT3", "DDIT4", "GADD45A", "PPP1R15A", "SARS",
          "SLC1A4", "TRIB3", "WARS", "YARS", "ABCF2", "MTHFR", "PON2", "LMO4", "CBX4")
ATF6 <- c("CALR", "CRELD2", "DNAJB11", "HERPUD1", "HSPA5", "HYOU1", "MANF", "MIS12", "PDIA4", "PDIA6", "SEL1L", "SLC39A14", "TMEM50B", "UGDH", "PLEKHA6", "STARD4")
XBP1 <- c("DERL2", "DNAJB9", "DNAJC10", "EDEM2", "HSPA13", "LMAN1", "OSTC", "PDIA5", "PLPP5", "SEC23B", "SEC31A", "SEC61A1", "SRP19", "SRPRB", "SSR1", "SSR3", "UFM1", "FICD", "MBNL2")


AT2 <- AddModuleScore(
  object = AT2,
  features = list(XBP1),
  name = "XBP1_score"
)

g <- VlnPlot(
  object = AT2,
  features = "XBP1_score1",   # AddModuleScore appends "1" to the name
  group.by = "SCT_snn_res.0.2",
  pt.size = 0.0)
g

ggsave(
  filename = "XBP1_Clusters.pdf",  
  plot = g,
  width = 6, height = 5
)

#Figure 3J: PERK Pathway Dotplot Clusters 
p <- DotPlot(
  object = AT2,
  features = PERK
) +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme_classic() +
  labs(
    x = "SCT_snn_res.0.2 cluster",
    y = "Gene"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
p

ggsave(
  filename = "PERKgenes_clusters_dotplot.pdf",
  plot = p,
  width = 8,
  height = 4
)
###############################################################################
#Figure 5 Analysis
##################################################################################
AT2 <- readRDS("mouse.rds")

#5A: UMAP by genotype
umap <- DimPlot(
  object = AT2,
  reduction = "umap",
  group.by = "genotype",
  pt.size = 0.5  
) +
  scale_color_manual(
    values = c(
      "ZZ" = "#F525E4",
      "WT" = "#0DD439",
      "MM" = "#423BFF"
    )
  ) +
  theme_classic() +     
  theme(
    plot.title = element_blank(),           
    axis.text = element_blank(),            
    axis.ticks = element_blank(),           
    legend.title = element_blank()          
  )
umap 

ggsave(
  filename = "mouse_UMAP_bygenotype.pdf",  
  plot = umap,
  width = 6, height = 5
)

#5B: UMAP by cell type 
umap <- DimPlot(
  object = AT2,
  reduction = "umap",
  group.by = "annot_supervised",
  pt.size = 0.5  
)  +
  theme_classic() +     
  theme(
    plot.title = element_blank(),           
    axis.text = element_blank(),            
    axis.ticks = element_blank(),           
    legend.title = element_blank()          
  )
umap

ggsave(
  filename = "Mouse_CellLineage_UMAP.pdf",  
  plot = umap,
  width = 8, height = 6
)

#5C Dotplot epithelial markers 
Idents(AT2) <- "annot_supervised"
lineages_to_plot <- c("AT1", "AT2", "basal", "ciliated", "secretory")
AT2_sub <- subset(AT2, idents = lineages_to_plot)

gene_list <- c(
  "Cav1", "Pdpn", "Hopx", "Sftpc", "Sftpb", "Abca3", "Lamp3", "Napsa","Scgb3a2", "Scgb3a1", "Scgb1a1", "Krt5", "Krt17","Trp63","Foxj1", "Dynlrb2","Ccdc153"
)
genes_to_plot <- gene_list[gene_list %in% rownames(AT2_sub)]

dotplot <- DotPlot(
  object = AT2_sub,
  features = genes_to_plot,
  group.by = "annot_supervised",
  assay = "SCT"
) +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme_classic() +
  RotatedAxis()

dotplot

ggsave(
  filename = "mouse_dotplot_lineages.pdf",  
  plot = dotplot,
  width = 8, height = 5)

#5D eYFP UMAP
umap <- FeaturePlot(
  object = AT2,
  feature = "eYFP",
  pt.size = 0.5  
)  +
  theme_classic() +     
  theme(
    plot.title = element_blank(),           
    axis.text = element_blank(),            
    axis.ticks = element_blank(),           
    legend.title = element_blank()          
  )
umap

ggsave(
  filename = "Mouse_eYFP_UMAP.pdf",  
  plot = umap,
  width = 6, height = 5
)

#5E Volcano Plot
AT2_sub <- subset(AT2, idents = "AT2")
Idents(AT2_sub) <- "genotype"

AT2_sub <- PrepSCTFindMarkers(AT2_sub)
DEG <- FindMarkers(
  AT2_sub,
  ident.1 = "ZZ",
  ident.2 = "MM",
  only.pos = FALSE,
  min.pct = 0.25,
  test.use = "wilcox"
)
write.csv(DEG, "Mouse_ZZvsMM_DEG.csv", row.names = TRUE)
sum(DEG$p_val_adj < 0.05 & DEG$avg_log2FC < -0.25,na.rm = TRUE) 

genes_ZZMM <- c("Egr1", "Ier3", "Junb", "Nfkbia", "Ier2", "Ier5", "Zfp36", "Btg2", "Hbegf",
                "Nr4a1", "Cdkn2a", "Atf3", "Jun", "Gadd45a", "Cebpd", "Dusp1", "Fos", "Rel", 
                "Phlda1", "Klf6", "Rhob", "Sat1", "Fosb")
genes_ZZWT <- c("eYFP", "Cdkn1a", "Hspa5", "Calr", "Hsp90b1", "Fosb", "Dnajb11", "Egr1", "Rel", "Ier3", "Ier5", "Nfkbia", "Irf1", "Hsbpb", "Herpud1", "Dnajb4", "Nfkb1", "Junb", "Nfkbiz")
genes_MMWT <- c("eYFP", "Tead1", "JAK2", "MAP4K3", "Jun", "Nfkbia", "Ier3", "Ier5", "Itga9", "Tomm6", "Nedd4", "Smad2", "Smad3", "Egr1", "Abca3", "Bmp1")

volcano <- EnhancedVolcano(
  DEG,
  lab = rownames(DEG),
  selectLab = genes_MMWT, 
  x = "avg_log2FC",
  y = "p_val_adj",
  pCutoff = 0.05,
  FCcutoff = 0.25,
  pointSize = 1.5,
  labSize = 5.0,
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  colAlpha = 0.8,
  legendPosition = 'right',
  legendLabSize = 8,
  legendIconSize = 3.0,
  title = 'ZZ vs MM DEGs',
  subtitle = NULL,
  caption = NULL,
  ylab = bquote(~-Log[10]~ 'padj'),
  xlab = bquote(~Log[2]~ 'FC'),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey30',
  max.overlaps = 10
) + 
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )
volcano

ggsave(
  filename = "Mouse_Volcano_MMvsWT.pdf",
  plot = volcano,
  width = 10, height = 6)

#5F Hallmark FGSEA
msig.H <- msigdbr(species = "Mus musculus", category = "H")
msig.H <- split(msig.H$gene_symbol, msig.H$gs_name)

gene_ranks <- DEG$avg_log2FC
names(gene_ranks) <- rownames(DEG)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

fgsea_res <- fgsea(
  pathways = msig.H,
  stats = gene_ranks,
  minSize = 5,
  maxSize = 500,
  scoreType = "std"
)
fgsea_res_df <- fgsea_res %>%
  mutate(leadingEdge = sapply(leadingEdge, paste, collapse = ";"))
write.csv(fgsea_res_df, "Mouse_MMvsWT_Hallmark.csv", row.names = FALSE)
fgsea_res$Significance <- factor(
  ifelse(fgsea_res$padj < 0.05, "padj < 0.05", "Not sig."),
  levels = c("padj < 0.05", "Not sig.")
)
fgsea_res$`-log(P-value)` <- -log10(fgsea_res$pval)
fgsea_res$pretty_pathway <- gsub("HALLMARK_", "", fgsea_res$pathway)
fgsea_res$pretty_pathway <- gsub("_", " ", fgsea_res$pretty_pathway)
fgsea_res_ord <- fgsea_res[order(abs(fgsea_res$NES), decreasing = TRUE), ]
fgseaRes_top <- head(fgsea_res_ord, 20)

p <- ggplot(
  fgseaRes_top,
  aes(
    x = reorder(str_wrap(str_sub(pretty_pathway, 1, 80), 40), NES),
    y = NES,
    fill = Significance
  )
) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  ggtitle("HALLMARK GSEA (MM vs WT)") +
  theme(
    axis.text = element_text(size = 7, face = "bold"),
    axis.title.y = element_blank()
  ) +
  scale_fill_manual(
    name = "Significance",
    values = c("padj < 0.05" = "blue", "Not sig." = "grey")
  ) +
  ylab("Normalized Enrichment Score (NES)")
p

ggsave(
  filename = "Mouse_Hallmark_MMvsWT.pdf",  
  plot = p,
  width = 6, height = 7
)

#5G: Venn Diagram Shared Genes 
genes1 <- read.csv("Penn_AATDvsCOPD_DEG.csv")
genes2 <- read.csv("LTRC_AATDvsCOPD_DEG.csv")
genes3 <- read.csv("iAT2_ZZvsMM_DEG.csv")
genes4 <- read.csv("Mouse_ZZvsMM_DEG.csv")
genes5 <- read.csv("Mouse_ZZvsWT_DEG.csv")

filter_deg <- function(df) {
  df %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC < -0.25) %>%
    dplyr::pull(X) %>%
    unique()
}

genes1_filtered <- filter_deg(genes1)
genes2_filtered <- filter_deg(genes2)
genes3_filtered <- filter_deg(genes3)
genes4_filtered <- filter_deg(genes4)
genes5_filtered <- filter_deg(genes5)

human_bulk_combined <- unique(c(
  genes1_filtered,
  genes2_filtered
))

mouse_to_human <- function(mouse_genes) {
  hg <- homologene(genes = mouse_genes, inTax = 10090, outTax = 9606)
  if(is.null(hg) || nrow(hg) == 0) {
    warning("No homologs found for these genes")
    return(character(0))
  }
  human_genes <- unique(as.character(hg$`9606`))
  unmapped <- setdiff(mouse_genes, hg$`10090`)
  if(length(unmapped) > 0) {
    message(length(unmapped), " mouse genes could not be mapped to human.")
  }
  return(human_genes)
}

genes4_human <- mouse_to_human(genes4_filtered)
manual_map <- c(
  "B230219D22Rik" = NA,    # no known human ortholog
  "Rab9"          = "RAB9A",
  "Eps8l2"        = "EPS8L2",
  "H2afy"         = "H2AFY",
  "Ralbp1"        = "RALBP1",
  "Klf3"          = "KLF3",
  "Trappc6b"      = "TRAPPC6B",
  "Atp6ap2"       = "ATP6AP2",
  "Smim6"         = "SMIM6",
  "Nsa2"          = "NSA2",
  "Mef2d"         = "MEF2D",
  "Psmd2"         = "PSMD2",
  "Cct5"          = "CCT5",
  "Glo1"          = "GLO1",
  "9530068E07Rik" = NA,
  "Pla2g1b"       = "PLA2G1B",
  "AU021092"      = NA,
  "Ppp1r2"        = "PPP1R2",
  "Psma6"         = "PSMA6",
  "Spry1"         = "SPRY1",
  "Cct2"          = "CCT2",
  "Wasf2"         = "WASF2",
  "Stard10"       = "STARD10",
  "Acaa1a"        = "ACAA1",
  "Matn4"         = "MATN4",
  "Minos1"        = "MINOS1",
  "Ppp2r2d"       = "PPP2R2D",
  "Marc2"         = "MARC2",
  "Gstt3"         = "GSTT3",
  "Ndel1"         = "NDEL1",
  "Ppig"          = "PPIG",
  "Arpc1b"        = "ARPC1B",
  "Epn3"          = "EPN3",
  "Stub1"         = "STUB1",
  "Spr"           = "SEPX1",
  "Pip4k2c"       = "PIP4K2C",
  "Rpn2"          = "RPN2",
  "Rtraf"         = "RTRA",
  "Zranb2"        = "ZRANB2",
  "Cald1"         = "CALD1",
  "2310030G06Rik" = NA,
  "Cct3"          = "CCT3",
  "Tsn"           = "TSN",
  "Arpc1a"        = "ARPC1A",
  "Vdac2"         = "VDAC2",
  "Fubp1"         = "FUBP1",
  "Cers2"         = "CERS2",
  "Kcnk1"         = "KCNK1",
  "Sec61a1"       = "SEC61A1",
  "Vdac1"         = "VDAC1",
  "Smap1"         = "SMAP1",
  "Ube2n"         = "UBE2N",
  "Foxa2"         = "FOXA2",
  "Krcc1"         = "KRCC1",
  "Cd47"          = "CD47",
  "Lrrc58"        = "LRRC58",
  "Dnajb1"        = "DNAJB1",
  "Slc48a1"       = "SLC48A1",
  "Rab6a"         = "RAB6A",
  "Rp9"           = "RP9",
  "Lbp"           = "LBP",
  "Rexo2"         = "REXO2",
  "Cbx4"          = "CBX4",
  "Chmp2a"        = "CHMP2A",
  "Suclg1"        = "SUCLG1",
  "Sumo3"         = "SUMO3",
  "Slc9a3r2"      = "SLC9A3R2",
  "Bmpr2"         = "BMPR2",
  "Rpn1"          = "RPN1",
  "Ostf1"         = "OSTF1",
  "Spty2d1"       = "SPTY2D1",
  "Pfdn2"         = "PFDN2",
  "Slc15a2"       = "SLC15A2",
  "Esd"           = "ESD",
  "Lgals3"        = "LGALS3",
  "Psmc6"         = "PSMC6",
  "Atp5o.1"       = "ATP5O",
  "Ccng1"         = "CCNG1",
  "Usp16"         = "USP16"
)
genes4_human <- genes4_human
for (gene in names(manual_map)) {
  if(!is.na(manual_map[gene])) {
    genes4_human <- c(genes4_human, manual_map[gene])
  }
}

genes5_human <- mouse_to_human(genes5_filtered)

manual_map2 <- c(
  "eYFP"          = NA,           # reporter, no human ortholog
  "Gt(ROSA)26Sor" = NA,           # mouse locus, no human ortholog
  "mt-Nd3"        = "MT-ND3",
  "mt-Atp8"       = "MT-ATP8",
  "Rps27"         = "RPS27",
  "Gm26561"       = NA,
  "Gm15564"       = NA,
  "Mt1"           = "MT1",
  "Acot1"         = "ACOT1",
  "Eloc"          = "ELOC",
  "Emsy"          = "EMSY",
  "Nsd3"          = "NSD3",
  "Fam208a"       = "FAM208A",
  "Sik1"          = "SIK1",
  "Hist1h1e"      = "HIST1H1E",
  "Gm20275"       = NA,
  "Cop1"          = "COP1",
  "Zfp871"        = NA,
  "Sult1a1"       = "SULT1A1",
  "Gm26917"       = NA,
  "Fam91a1"       = "FAM91A1",
  "Nectin3"       = "NECTIN3",
  "Selenot"       = "SELENOT"
)
genes5_human <- genes5_human
for (gene in names(manual_map2)) {
  if(!is.na(manual_map2[gene])) {
    genes5_human <- c(genes5_human, manual_map2[gene])
  }
}

mouse_combined_human <- unique(c(
  genes4_human,
  genes5_human
))

gene_lists <- list(
  Human_bulk_DEGs = human_bulk_combined,
  ZZ_vs_MM_iAT2  = genes3_filtered,
  Mouse_DEGs     = mouse_combined_human
)

venn.plot <- venn.diagram(
  x = gene_lists,
  filename = NULL,              
  fill = c("purple", "blue", "lightgreen")[1:length(gene_lists)],
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.0,
  cat.pos = 0,
  margin = 0.1,
  main = "Gene Set Overlap"
)

grid.newpage()
grid.draw(venn.plot)

pdf("gene_venn.pdf", width = 6, height = 6)
grid.draw(venn.plot)
dev.off()

long_df <- enframe(gene_lists, name = "Set", value = "Gene") %>%
  unnest(Gene)

gene_matrix <- long_df %>%
  mutate(Present = 1) %>%
  pivot_wider(names_from = Set, values_from = Present, values_fill = 0)

gene_matrix <- gene_matrix %>%
  unite("Intersection", names(gene_lists), sep = "_", remove = FALSE)

disjoint_genes <- gene_matrix %>%
  rowwise() %>%
  mutate(Sets_present = sum(c_across(names(gene_lists)))) %>%
  filter(Sets_present >= 1) %>%
  ungroup()

intersection_lists <- disjoint_genes %>%
  group_by(Intersection) %>%
  summarise(Genes = paste(Gene, collapse = ";"))

write.csv(intersection_lists, "venn_disjoint_intersections_down.csv", row.names = FALSE)

#Figure 5I: Dotplot Hallmark Pathways across datasets
files <- c( "Penn_Hallmark.csv", "LTRC_Hallmark.csv", "iAT2_Hallmark.csv", "Mouse_ZZvsWT_Hallmark.csv", "Mouse_ZZvsMM_Hallmark.csv" ) 
dataset_names <- c("Penn", "LTRC", "iAT2", "Mouse_ZZvsWT", "Mouse_ZZvsMM") 
my_pathways <- c( "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_APOPTOSIS","HALLMARK_P53_PATHWAY" )

fgsea_long <- purrr::map2_df(files, dataset_names, ~ {
  read.csv(.x) %>%
    tibble::as_tibble() %>%
    mutate(dataset = .y)
})

fgsea_manual <- fgsea_manual %>%
  filter(pathway %in% my_pathways) %>%
  mutate(
    pathway = factor(pathway, levels = rev(my_pathways)),  
    dataset = factor(dataset, levels = c("Penn", "LTRC", "iAT2", "Mouse_ZZvsMM", "Mouse_ZZvsWT")), # x-axis order
    dot_size = -log10(padj)
  )

g <- ggplot(fgsea_manual, aes(x = dataset, y = pathway)) +
  geom_point(aes(fill = NES, size = dot_size), shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "NES"
  ) +
  scale_size(range = c(3, 8), name = "-log10(padj)") +
  theme_minimal(base_size = 12) +
  labs(
    x = NULL,
    y = NULL,
    title = "FGSEA Hallmark DotPlot Across Datasets"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(color = "gray90")
  )
g

ggsave(
  filename = "Hallmark_all.pdf",  
  plot = g,
  width = 8, height = 10
)

