#Load libraries
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(EnhancedVolcano)
  library(msigdbr)
  library(fgsea)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(dorothea)
  library(viper)
  library(pheatmap)
  library(VennDiagram)
  library(enrichR)
  library(homologene)

######################################
#Figure 2
AT2 <- readRDS("iAT2.rds")

# Figure 2B: UMAP by genotype
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

#Figure 2C: Louvain Clusters
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

#Figure 2D: AT2 & non-AT2 markers 
Idents(AT2) <- "SCT_snn_res.0.2"
genes <- c("NKX2-1","SFTPB","SLC34A2","ABCA3", #AT2
           "CDX2", #intestinal
           "AFP",#liver
           "TP63", #basal
           "SCGB3A2", #secretory
           "FOXJ1" #ciliated)

p <- DotPlot(
  object = AT2,
  features = genes
) +
  scale_color_gradient(
    low = "lightgrey",
    high = "blue"
  ) +
  scale_size(range = c(.25, 10)) +   
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) +
  labs(
    x = "Gene",
    y = "Cluster"
  )
p

ggsave(
  filename = "iAT2_clusters.pdf",
  plot = p,
  width = 8,
  height = 4
)

#Figure 2E: Heatmap of DEGs by cluster
Idents(AT2) <- "SCT_snn_res.0.2"
DEG <- FindAllMarkers(
  AT2,
  only.pos = TRUE,
  min.pct = 0.10,
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

#Figure 2F & L: Cluster Defining Genes UMAP 
umap <- FeaturePlot(AT2, feature = "HSPA6") #cluster 4
umap <- FeaturePlot(AT2, feature = "IFI44L") #cluster 1
umap <- FeaturePlot(AT2, feature = "TSPO") #cluster 0
umap <- FeaturePlot(AT2, feature = "TOP2A") #clusters 2/3
umap <- FeaturePlot(AT2, feature = "KRT17") #ABI

ggsave(
  filename = "KRT17_UMAP.pdf",  
  plot = umap,
  width = 6, height = 5
)

#Figure 2G: Hallmark Pathways Defining Each Cluster
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

#Figure 2H: Heatmap directly comparing top DEGs clusters 4 & 1
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

#Figure 2I: ER stress specific GO pathways cluster 4 vs 1
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

#Figure 2J & S2H/I: UPR arms module scores 
PERK <- c("ASNS", "ATF3", "ATF4", "ATF6", "CARS", "CEBPB", "DDIT3", "DDIT4", "GADD45A", "PPP1R15A", "SARS",
          "SLC1A4", "TRIB3", "WARS", "YARS", "ABCF2", "MTHFR", "PON2", "LMO4", "CBX4")
ATF6 <- c("CALR", "CRELD2", "DNAJB11", "HERPUD1", "HSPA5", "HYOU1", "MANF", "MIS12", "PDIA4", "PDIA6", "SEL1L", "SLC39A14", "TMEM50B", "UGDH", "PLEKHA6", "STARD4")
XBP1 <- c("DERL2", "DNAJB9", "DNAJC10", "EDEM2", "HSPA13", "LMAN1", "OSTC", "PDIA5", "PLPP5", "SEC23B", "SEC31A", "SEC61A1", "SRP19", "SRPRB", "SSR1", "SSR3", "UFM1", "FICD", "MBNL2")


AT2 <- AddModuleScore(
  object = AT2,
  features = list(PERK),
  name = "PERK_score"
)

g <- VlnPlot(
  object = AT2,
  features = "PERK_score1",   # AddModuleScore appends "1" to the name
  group.by = "SCT_snn_res.0.2",
  pt.size = 0.0)
g

ggsave(
  filename = "PERK_Clusters.pdf",  
  plot = g,
  width = 6, height = 5
)

#Figure 2K: PERK Pathway Dotplot Clusters 
p <- DotPlot(
  object = AT2,
  features = PERK
) +
  scale_color_gradient(
    low = "lightgrey",
    high = "blue"
  ) +
  scale_size(range = c(.25, 10)) +   
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) +
  labs(
    x = "Gene",
    y = "Cluster"
  )
p

ggsave(
  filename = "PERKgenes_clusters_dotplot.pdf",
  plot = p,
  width = 10,
  height = 5
)

#Figure 2M: ABI module score
ABI <- c("CDH1","CDH2","SPINK1","MMP7","PTGS2","CDKN2A","CDKN2B","HMGA2","EPCAM","VIM","FN1","COL1A1","TNC","VCAN","PCP4","CUX2","PRSS2","CPA6","CTSE","MDK","GDF15","SLCO2A1","EPHB2","ITGB8","ITGAV","ITGB6","TGFBI","KCNN4","KCNQ5","KCNS3","CDKN1A","CCND1","CCND2","MDM2","OCIAD2","PTCHD4")

AT2 <- AddModuleScore(
  object = AT2,
  features = list(ABI),
  name = "ABI_score"
)

p <- VlnPlot(
  object = AT2,
  features = "ABI_score1",
  group.by = "SCT_snn_res.0.2",
  pt.size = 0
) +
  theme_classic() +
  ylab("ABI module score") +
  xlab("Cluster")
p

ggsave(
  filename = "iAT2_ABI.pdf",  
  plot = p,
  width = 6, height = 5
)

#Figure 2N: ABI Dotplot
Idents(AT2) <- "SCT_snn_res.0.2"
top20 <- c("CDH1","CDH2","SPINK1","MMP7","PTGS2","CDKN2A","CDKN2B","HMGA2","EPCAM","VIM","FN1","COL1A1","TNC","VCAN","PCP4","CUX2","PRSS2","CPA6","CTSE","MDK")
           
p <- DotPlot(
  object = AT2,
  features = top20
) +
  scale_color_gradient(
    low = "lightgrey",
    high = "blue"
  ) +
  scale_size(range = c(.25, 8)) +   
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) +
  labs(
    x = "Gene",
    y = "Cluster"
  )
p

ggsave(
  filename = "iAT2_ABI.pdf",
  plot = p,
  width = 8,
  height = 4
)

#Supplement 2A
Idents(AT2) <- "orig.ident"
DEG <- FindMarkers(
  AT2,
  ident.1 = "PiZZ_Air",
  ident.2 = "PiMM_Air",
  only.pos = FALSE,
  min.pct = 0.10,
  test.use = "wilcox"
)

write.csv(DEG, "iAT2_ZZvsMM_DEG.csv", row.names = TRUE)
sum(DEG$p_val_adj < 0.05 & DEG$avg_log2FC > 0.25, na.rm = TRUE) # < -.25 for COPD

geneslabel <- c("COL4A2", 'HSPA1A', "OAS2", 'IFITM1', 'IFI44L', "DDIT3", "ATF3", "TNF", "DNAJB1", "GADD45B")
  
volcano <- EnhancedVolcano(
  DEG,
  lab = rownames(DEG),
  x = "avg_log2FC",
  y = "p_val_adj",
  selectLab = geneslabel,
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
  max.overlaps = 10, 
) + 
  theme(
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

volcano

ggsave(
  filename = "iAT2_Volcano_ZZvsMM.pdf",
  plot = volcano,
  width = 10, height = 6)

#Supplement 2B
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

write.csv(fgsea_res_df, "iAT2_Hallmark.csv", row.names = FALSE)

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
  filename = "iAT2_Hallmark_ZZvsMM.pdf",  
  plot = p,
  width = 6, height = 7
)

#################################################
#Figure 4

AT2 <- readRDS("mouse.rds")

#4B: UMAP by genotype
emb <- Embeddings(AT2, "umap") %>% as.data.frame()
emb$cell <- rownames(emb)
meta <- AT2@meta.data %>% as.data.frame()
meta$cell <- rownames(meta)
df <- left_join(emb, meta, by = "cell")
set.seed(123)
df_shuffled <- df %>% slice_sample(prop = 1)

umap <-  ggplot(df_shuffled, aes(x = umap_1, y = umap_2, color = genotype)) +
  geom_point(size = 0.5, alpha = 0.85) +
  scale_color_manual(values = c(
    "ZZ" = "#F525E4",
    "WT" = "#0DD439",
    "MM" = "#423BFF"
  )) +
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

#4C: UMAP by cell type 
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

#4D: Dotplot epithelial markers 
Idents(AT2) <- "annot_supervised"
lineages_to_plot <- c("AT1", "AT2", "basal", "ciliated", "secretory")
AT2_sub <- subset(AT2, idents = lineages_to_plot)

gene_list <- c(
  "Cav1", "Pdpn", "Hopx", #AT1
  "Sftpc", "Sftpb", "Abca3", "Lamp3", "Napsa", #AT2
  "Scgb3a2", "Scgb3a1", "Scgb1a1", #secretory
  "Krt5", "Krt17", #ABI
  "Trp63", #Basal
  "Foxj1", "Dynlrb2","Ccdc153" #ciliated
)
genes_to_plot <- gene_list[gene_list %in% rownames(AT2_sub)]

p <- DotPlot(
  object = AT2_sub,
  features = genes_to_plot
) +
  scale_color_gradient(
    low = "lightgrey",
    high = "blue"
  ) +
  scale_size(range = c(.25, 8)) +   
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) +
  labs(
    x = "Gene",
    y = "Cluster"
  )
p

ggsave(
  filename = "mouse_dotplot_lineages.pdf",  
  plot = p,
  width = 8, height = 5)

#4E: eYFP UMAP
eyfp_df <- FetchData(AT2, vars = "eYFP") %>% as.data.frame()
eyfp_df$cell <- rownames(eyfp_df)
colnames(eyfp_df)[colnames(eyfp_df) == "eYFP"] <- "eYFP_expr"
df <- left_join(df, eyfp_df, by = "cell")
set.seed(123)
cells_shuffled <- sample(colnames(AT2))
eYFP <-FeaturePlot(AT2, features = "eYFP", cells = cells_shuffled)

ggsave(
  filename = "Mouse_eYFP_UMAP.pdf",  
  plot = eYFP,
  width = 6, height = 5
)

#4F & S4B/C: Volcano Plot DEGs 
Idents(AT2) <- "annot_supervised"
AT2_sub <- subset(AT2, idents = "AT2")
Idents(AT2_sub) <- "genotype"

AT2_sub <- PrepSCTFindMarkers(AT2_sub)
DEG <- FindMarkers(
  AT2_sub,
  ident.1 = "ZZ",
  ident.2 = "WT",
  only.pos = FALSE,
  min.pct = 0.10,
  test.use = "wilcox"
)
write.csv(DEG, "Mouse_ZZvsWT_DEG.csv", row.names = TRUE)
sum(DEG$p_val_adj < 0.05 & DEG$avg_log2FC > 0.25,na.rm = TRUE) 

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

#4G & S4D/E: Hallmark FGSEA
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

#######################################################
#Figure 5
AT2_1 <- readRDS("LTRC_AT2only.rds")

# Figure 5B: LTRC Volcano plot
Idents(AT2_1) <- "disease_state"

DEG <- FindMarkers(
  AT2_1,
  ident.1 = "AATD",
  ident.2 = "COPD",
  only.pos = FALSE,
  min.pct = 0.1,
  test.use = "wilcox", 
  recorrect_umi = FALSE
)
write.csv(DEG, "LTRC_AATDvsCOPD_DEG.csv", row.names = TRUE)

sum(DEG$p_val_adj < 0.05 & DEG$avg_log2FC < -0.25, na.rm = TRUE) # < -.25 for COPD

genes_lable <- c("NFKBIZ","PELI1","CXCL2","SOCS3","LIFR","AREG","IL34", "FOS","FOSB","JUN","JUNB","JUND","ATF3","EGR1","NR4A1","NR4A2","NR4A3","KLF6","BTG2","CSRNP1","JDP2", "XBP1","HSPA5","EDEM1","TXNDC11","MANBA","ATF3","GADD45B","TP53INP1","SESN3","DEPTOR","ULK2", "ZFAND5")

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
  legendLabSize = 3,
  legendIconSize = 3.0,
  title = NULL,
  subtitle = NULL,
  caption = NULL,
  ylab = bquote(~-Log[10]~ 'padj'),
  xlab = bquote(~Log[2]~ 'FC'),
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = 'grey30',
  max.overlaps = 35, 
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

#Figure 5C: LTRC Hallmark GSEA 
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

write.csv(fgsea_res_df, "LTRC_AATDvsCOPD_Hallmark.csv", row.names = FALSE)

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

#Figure 5D: PROGENy 
pathway_mat <- progeny::getModel(organism = "Human", top = 500)

expr_mat <- GetAssayData(AT2_1, layer = "counts", assay = "SCT")

common_genes <- intersect(rownames(expr_mat), rownames(pathway_mat))
cat("Matched genes:", length(common_genes), "out of", 
    nrow(pathway_mat), "pathways genes\n")

pathway_scores <- t(pathway_mat[common_genes, ]) %*% expr_mat[common_genes, ]

pathway_scores <- t(scale(t(pathway_scores)))

AT2_1@meta.data <- cbind(
  AT2_1@meta.data,
  as.data.frame(t(pathway_scores))
)

head(colnames(AT2_1@meta.data))

pathway_cols <- rownames(pathway_scores)
pathway_means <- AT2_1@meta.data %>%
  group_by(disease_state) %>%
  summarise(across(all_of(pathway_cols), mean, na.rm = TRUE), .groups = "drop")

pathway_heatmap_mat <- as.matrix(pathway_means[, -1])
rownames(pathway_heatmap_mat) <- pathway_means$disease_state

pathway_aatd_mean <- pathway_means %>%
  filter(disease_state == "AATD") %>%
  pivot_longer(cols = -disease_state, names_to = "pathway", values_to = "score") %>%
  arrange(desc(score))

pathway_heatmap_mat_ordered <- pathway_heatmap_mat[, pathway_aatd_mean$pathway]

p <- pheatmap::pheatmap(
  t(pathway_heatmap_mat_ordered),
  main = "Pathway Scores: AATD vs COPD (ordered by AATD enrichment)",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  display_numbers = FALSE,  
  fontsize_number = 10,
  cluster_cols = FALSE,     
  cluster_rows = FALSE,  
  filename = "LTRC_AT2_Progeny_AATDvsCOPD.pdf", 
  width = 6, 
  height = 8.5,
  cellwidth = 60,
  cellheight = 15
)

#S5A: Volcano plot UPenn
AT2_2 <- readRDS("AATD.Rds")
Idents(AT2_2) <- "annot3"
AT2_subset <- subset(AT2_2, idents = "AT2")
Idents(AT2_subset) <- "annot2"

DEG <- FindMarkers(
  AT2_subset,
  ident.1 = "AATD",
  ident.2 = "COPD",
  only.pos = FALSE,
  min.pct = 0.10,
  test.use = "wilcox", 
  recorrect_umi = FALSE
)
write.csv(DEG, "Penn_AATDvsCOPD_DEG.csv", row.names = TRUE)

sum(DEG$p_val_adj < 0.05 & DEG$avg_log2FC > 0.25,na.rm = TRUE) # < -.25 for COPD

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

#Figure 5E: Venn diagram & EnrichR shared genes LTRC & UPenn
  deg_ltrc <- read.csv("LTRC_AATDvsCOPD_DEG.csv", row.names = 1) 
deg_penn <- read.csv("Penn_AATDvsCOPD_DEG.csv", row.names = 1) 

genes_ltrc <- rownames( deg_ltrc[deg_ltrc$avg_log2FC > 0.25 & deg_ltrc$p_val_adj < 0.05, ] ) 
genes_penn <- rownames( deg_penn[deg_penn$avg_log2FC > 0.25 & deg_penn$p_val_adj < 0.05, ] ) 
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
head(sig[, c("Term", "Adjusted.P.value", "Odds.Ratio", "Genes")], 25)

##################################################################
#Figure 6

#Figure 6A/B & S6A-C: Venn Diagram & GO Pathways
genes1 <- read.csv("Penn_AATDvsCOPD_DEG.csv")
genes2 <- read.csv("LTRC_AATDvsCOPD_DEG.csv")
genes3 <- read.csv("iAT2_ZZvsMM_DEG.csv")
genes4 <- read.csv("Mouse_ZZvsMM_DEG.csv")
genes5 <- read.csv("Mouse_ZZvsWT_DEG.csv")

filter_deg <- function(df) {
  df %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0.25) %>%
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

#Figure 6C: Hallmark dotplot
files <- c( "Penn_AATDvsCOPD_Hallmark.csv", "LTRC_AATDvsCOPD_Hallmark.csv", "iAT2_Hallmark.csv", "Mouse_ZZvsWT_Hallmark.csv", "Mouse_ZZvsMM_Hallmark.csv" ) 
dataset_names <- c("Penn", "LTRC", "iAT2", "Mouse_ZZvsWT", "Mouse_ZZvsMM") 
my_pathways <- c( "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_ESTROGEN_RESPONSE_EARLY", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_APOPTOSIS","HALLMARK_P53_PATHWAY" )

fgsea_long <- purrr::map2_df(files, dataset_names, ~ {
  read.csv(.x) %>%
    tibble::as_tibble() %>%
    mutate(dataset = .y)
})

fgsea_manual <- fgsea_long %>%
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

#Figure 6D & S6D: NFkB Module score 
NFkB <- c("NFKB1", "NFKB2", "NFKBIA", "RELA", "TRAF2", "CCL2", "ICAM1", "SERPINA3", "TNFAIP3", "CCL5")
NFkB_mouse <- c("Nfkb1", "Nfkb2", "Nfkbia", "Rela", "Traf2", "Ccl2", "Icam1", "Serpina3n", "Tnfaip3", "Ccl5")

#change based on sample: AT2 = iAT2, AT2_sub = mouse, AT2_1 = LTRC, AT2_subset = UPenn
Idents(AT2_1) <- "disease_state" 
AT2_1 <- subset(AT2_1, idents = c("AATD", "COPD"))
AT2_1 <- AddModuleScore(
  object = AT2_1,
  features = list(NFkB),
  name = "NFkB_score"
)

p <- VlnPlot(
  object = AT2_1,
  features = "NFkB_score1",
  group.by = "disease_state",
  pt.size = 0
) +
  scale_fill_manual(
    values = c(
      "AATD" = "#F525E4",
      "COPD" = "#423BFF",
      "WT" = "green"
    )
  ) +
  theme_classic() +
  ylab("NFkB module score") +
  xlab("Condition")

p

ggsave(
  filename = "NFkB_LTRC.pdf",  
  plot = p,
  width = 5, height = 6
)
#Figure 6E & S6E: Inflammatory AT2 module score 
AT2i <- c("ITGA2", "TNC", "CXCL8", "TNIP3", "RASGRP1", "SGPP2", "CSF3", "TNFAIP3", "NABP1", "DDR2", "LAMB3", "BATF", "STAT4", "CD83", "EHD1", "ICAM4", "SLC6A14", "CNKSR3", "RND1", "ICAM5", "ABCA13", "TXNRD1", "IPCEF1", "S100A10", "PFKP",
          "HIVEP2", "ELOVL7", "TNFRSF21", "TNFSF14", "NFKB1", "TNIP1", "DMBT1", "IRAK2", "LYN", "BIRC3", "PELI1", "SLC5A8", "ANXA2", "TNFRSF10D", "NAMPT", "UAP1", "ICAM1", "MTHFD2L", "AGPAT5", "HIF1A", "DGKH", "SAMD4A", "EID3", "SLC7A1")
AT2i_mouse <- c("Itga2","Tnc","Cxcl15","Tnip3","Rasgrp1","Sgpp2","Csf3","Tnfaip3","Nabp1","Ddr2","Lamb3","Batf","Stat4","Cd83","Ehd1","Icam4","Slc6a14","Cnksr3","Rnd1","Icam5","Abca13","Txnrd1","Ipcef1","S100a10","Pfkp","Hivep2","Elovl7","Tnfrsf21","Tnfsf14","Nfkb1","Tnip1","Dmbt1","Irak2","Lyn","Birc3","Peli1","Slc5a8","Anxa2","Tnfrsf10d","Nampt","Uap1","Icam1","Mthfd2l","Agpat5","Hif1a","Dgkh","Samd4a","Eid3","Slc7a1")

#change based on sample: AT2 = iAT2, AT2_sub = mouse, AT2_1 = LTRC, AT2_subset = UPenn
AT2_1 <- AddModuleScore(
  object = AT2_1,
  features = list(AT2i),
  name = "AT2i_score"
)

p <- VlnPlot(
  object = AT2_1,
  features = "AT2i_score1",
  group.by = "disease_state",
  pt.size = 0
) +
  scale_fill_manual(
    values = c(
      "AATD" = "#F525E4",
      "COPD" = "#423BFF",
      "Control" = "green"
    )
  ) +
  theme_classic() +
  ylab("AT2i module score") +
  xlab("Condition")
p

ggsave(
  filename = "AT2i_iAT2.pdf",  
  plot = p,
  width = 5, height = 6
)
           
#Figure 6F: Inflammatory AT2 genes dotplot 
AT2i <- c("SGPP2", "NFKB1", "SOD2", "HIF1A", "CSF3", "IL4R", "TNFAIP3", "CXCL8", "CXCL1", "CCL20", "CXCL3")

#change based on sample: AT2 = iAT2, AT2_1 = LTRC, AT2_subset = UPenn
Idents(AT2_1) <- "disease_state"
p <- DotPlot(
  object = AT2_1,
  features = AT2i
) +
  scale_color_gradient(
    low = "lightgrey",
    high = "blue"
  ) +
  scale_size(range = c(.25, 8)) +   
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  
  ) +
  labs(
    x = "Gene",
    y = "Cluster"
  )
p

ggsave(
  filename = "AT2i_LTRC_Dotplot.pdf",
  plot = p,
  width = 8,
  height = 4
)
