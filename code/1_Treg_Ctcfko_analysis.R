rm(list = ls())
gc()

library(reticulate)
use_condaenv("scrna", required = TRUE)

library(tidyverse)
library(Seurat)
library(cowplot)
library(scop)
library(ggsci)
library(colorspace)
library(RColorBrewer)
library(patchwork)
library(ggsignif)
source("/home/data/qiuli/Bioinformatics/code2source/SCP_adjustment_for_v5.R", echo=TRUE)

alldata <- readRDS(file = "results/3_anno.rds")
CellDimPlot(alldata, group.by = "main_anno", reduction = "umap.scanorama", label = T)
CellDimPlot(alldata, group.by = "RNA_snn_res.0.8", reduction = "umap.scanorama", label = T)
ggsave(filename = "results/4_funtional_analysis/0_UMAP_main_anno.pdf", height = 5, width = 5)
alldata$main_anno <- factor(alldata$main_anno)
CellStatPlot(alldata, stat.by = "main_anno", group.by = "group", plot_type = "trend", label = T)
ggsave(filename = "results/4_funtional_analysis/0_main_anno_proportion.pdf", width = 3.5, height = 5)

#### CD8T ----
rm(list = ls()[-which(ls() %in% c("alldata"))])
gc()
CD8T <- subset(alldata, subset = main_anno == "CD8T")
CellDimPlot(CD8T, group.by = "fine_anno", reduction = "umap.scanorama", label = T)
CellDimPlot(CD8T, group.by = "RNA_snn_res.0.8", reduction = "umap.scanorama", label = T)
selected_cells <- Embeddings(CD8T, reduction = "umap.scanorama") %>% as.data.frame() %>% 
  filter(umapscanorama_2 >= -3 & umapscanorama_2 <= 6) %>% rownames()
CD8T <- subset(CD8T, cells = selected_cells)
DimPlot(CD8T, group.by = "RNA_snn_res.0.8", reduction = "umap.scanorama", label = T)
Idents(CD8T) <- "RNA_snn_res.0.8"
CD8T <- FindSubCluster(CD8T, cluster = 10, graph.name = "RNA_snn", subcluster.name = "celltype1", resolution = 0.3)
DimPlot(CD8T, group.by = "celltype1", reduction = "umap.scanorama", label = T)
CD8T$celltype1 <- factor(CD8T$celltype1)
levels(CD8T$celltype1) <- c("C2", "C2", "C1", "C3", "C6", "C4", "C3", "C4", "C5", "C3")
CD8T$celltype1 <- factor(CD8T$celltype1, levels = paste0("C", 1:6))
Idents(CD8T) <- "celltype1"
CD8_markers <- FindAllMarkers(CD8T)
write_csv(CD8_markers, file = "results/4_funtional_analysis/1_CD8T_markers.csv")

colors <- brewer.pal(12, "Paired")
CellDimPlot(CD8T, group.by = "celltype1", reduction = "umap.scanorama", label = T, palcolor = colors[6:1])
ggsave(filename = "results/4_funtional_analysis/1_CD8T_cluster.pdf", height = 5, width = 5)
CellDimPlot(CD8T, group.by = "celltype1", reduction = "umap.scanorama", split.by = "group", label = T, palcolor = colors[6:1])
ggsave(filename = "results/4_funtional_analysis/1_CD8T_cluster_by_group.pdf", height = 5, width = 10)
FeatureDimPlot_v5(CD8T, features = c("Tcf7", "Havcr2", "Sell", "Gzmb", "Cxcr6", "Cx3cr1"), reduction = "umap.scanorama")
ggsave(filename = "results/4_funtional_analysis/1_CD8T_gene_umap.pdf", height = 5, width = 10)
FeatureDimPlot_v5(CD8T, features = c("Cd44", "Sell", "Cxcr6", "Cx3cr1", "Fasl", "Gzmk"), reduction = "umap.scanorama")
ggsave(filename = "results/4_funtional_analysis/1_CD8T_gene_umap2.pdf", height = 5, width = 10)
FeatureDimPlot_v5(CD8T, features = c("Prf1", "Ifng"), reduction = "umap.scanorama")
ggsave(filename = "results/4_funtional_analysis/1_CD8T_gene_umap3.pdf", height = 3, width = 6)
FeatureDimPlot_v5(CD8T, features = c("Prf1", "Pdcd1", "Lag3", "Tigit", "Tox"), reduction = "umap.scanorama")
FeatureStatPlot_v5(CD8T, stat.by = c("Prf1", "Gzmb", "Pdcd1", "Lag3", "Tigit", "Tox"), group.by = "celltype1",
                   bg.by = "celltype1", stack = T)
Idents(CD8T) <- "celltype1"
StackedVlnPlot(obj = CD8T, features = c("Prf1", "Gzmb", "Pdcd1", "Lag3", "Tigit", "Tox"))
VlnPlot(CD8T, c("Prf1", "Gzmb", "Pdcd1", "Lag3", "Tigit", "Tox"), stack = TRUE, sort = FALSE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")

# https://www.nature.com/articles/s41467-023-35948-9#Bib1
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7341730/#R53
# https://pubmed.ncbi.nlm.nih.gov/29716982/
# https://www.science.org/doi/10.1126/science.aaz8658
# https://www.nature.com/articles/s41590-023-01685-w
FeatureDimPlot_v5(CD8T, features = c("Hspa1a", "Hspa1b", "Hsph1", "Dnajb1", "Klf2", "Foxp1"), reduction = "umap.scanorama", split.by = "group")
ggsave(filename = "results/4_funtional_analysis/1_CD8T_gene_umap_hsp.pdf", height = 8, width = 12)

### 1. heatmap
CD8T <- RunDEtest(CD8T, group_by = "celltype1")
markers <- CD8T@tools[["DEtest_celltype1"]][["AllMarkers_wilcox"]]
de_filter <- markers %>% 
  mutate(group1 = factor(group1, levels = paste0("C", 1:6))) %>% 
  filter(p_val_adj < 0.05 & avg_log2FC > 1 & pct.1 >= 0.3) %>% 
  group_by(group1) %>% 
  arrange(desc(avg_log2FC)) %>%   
  slice_head(n = 60) %>% 
  ungroup()
expr <- AverageExpression(CD8T, assays = "RNA", features = de_filter$gene, group.by = "celltype1")$RNA %>% as.matrix()
expr <- expr[de_filter$gene, paste0("C", 1:6)]
expr <- log2(expr + 1)
expr <- apply(expr, 1, scale) %>% t()
colnames(expr) <- paste0("C", 1:6)

library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(c(-2,-1.5,0,1.5,2), c("#1E4782", "#66A1CC", "white", "#FB9A99", "#B42D34"))
column_ha <- HeatmapAnnotation(Cluster = colnames(expr), 
                               col = list(Cluster = c("C1" = colors[6],"C2" = colors[5],"C3" = colors[4],
                                                      "C4" = colors[3],"C5" = colors[2],"C6" = colors[1])),
                               show_annotation_name = F)
pdf(file = "results/4_funtional_analysis/1_CD8T_heatmap.pdf", height = 8, width = 4)
Heatmap(expr, col = col_fun, top_annotation = column_ha, cluster_rows = F, cluster_columns = F, 
        column_split = paste0("C", 1:6), row_split = c(rep(1,60), rep(2,60), rep(3,38), rep(4,12), rep(5,33), rep(6,60)),
        show_row_names = F, show_column_names = T, column_names_side = "top",
        column_names_rot = 0)
dev.off()

### 1. Dotplot
genes2plot <- c("Cd44", "Sell", "Tcf7", "Slamf6", "Il7r", "Cxcr6", "Cx3cr1", "Gzmb", "Icos", "Prf1", "Pdcd1", "Havcr2", "Lag3", "Tigit", "Tox", "Stmn1", "Cdk1")
DotExpr <- AverageExpression(CD8T, assays = "RNA", features = genes2plot, group.by = "celltype1")$RNA %>% as.matrix()
DotExpr <- log2(DotExpr[genes2plot,] + 1)
Idents(CD8T) <- "celltype1"
markers1 <- FindAllMarkers(CD8T, features = genes2plot, logfc.threshold = 0, min.pct = 0)
DotData <- as.data.frame(DotExpr) %>% rownames_to_column("gene") %>% 
  gather(key = "cluster", value = "expression", -gene) %>% 
  left_join(., select(markers1, c(gene, cluster, pct.1)), by = c("gene", "cluster"))
DotData$pct.1 <- ifelse(is.na(DotData$pct.1), 0, DotData$pct.1)
DotData$gene <- factor(DotData$gene, levels = rev(genes2plot))

theme_niwot <- function(){
  theme_bw() +
    theme(
      text = element_text(family = "sans"),
      axis.line.x = element_line(color = "black"),
      axis.line.y = element_line(color = "black"),
      axis.text.x = element_text(family = "sans", size = 12, face = "plain", color = "black"),
      axis.text.y = element_text(family = "sans", size = 12, face = "plain", color = "black"),
      panel.border = element_blank(),
      axis.title.x = element_text(margin = margin(t = 10), size = 14,
                                  family = "sans", color = "black"),
      axis.title.y = element_text(margin = margin(t = 10), size = 14,
                                  family = "sans", color = "black"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      plot.margin = unit(c(1 ,1, 1, 1), units = "cm"),
      legend.text = element_text(size = 12, family = "sans"),
      legend.title = element_blank(),
      legend.key = element_blank(),
      panel.background = element_rect(fill = "white"),
      legend.background = element_rect(color = "black", fill = "transparent", linetype = "blank")
    )
}

ggplot(DotData, aes(x = cluster, y = gene)) +
  geom_point(aes(size = pct.1, fill = expression), color = "black", shape = 21, stroke = 0.5) +
  scale_fill_gradient2(low = "#F6F6F6", mid = "#FF888E", high = "#A01021", midpoint = 2.5) +
  theme_niwot() +
  labs(x = "", y = "")
ggsave(filename = "results/4_funtional_analysis/1_CD8T_dotplot.pdf", height = 4.5, width = 4.5)

### 2. Differential expression analysis
C3vsC4 <- FindMarkers(CD8T, ident.1 = "C3", ident.2 = "C4", logfc.threshold = 0) %>% rownames_to_column("gene")
write_csv(C3vsC4, file = "results/4_funtional_analysis/2_C3vsC4.csv")

Tsm <- subset(CD8T, subset = celltype1 %in% c("C2"))
Idents(Tsm) <- "group"
TsmCkovsCtr <- FindMarkers(Tsm, ident.1 = "cko", ident.2 = "ctr") %>% rownames_to_column("gene")
write_csv(TsmCkovsCtr, file = "results/4_funtional_analysis/2_TsmCkovsCtr.csv")

library(readxl)
C2 <- read_xlsx(path = "results/4_funtional_analysis/2_GSEA_C2_CP_all.xlsx")
C5 <- read_xlsx(path = "results/4_funtional_analysis/2_GSEA_C5_GO_all.xlsx")
GSEA_res <- rbind(C2, C5) %>% filter(p.adjust < 0.05)
GSEA_dot <- GSEA_res %>% 
  filter(ID %in% c("REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM", "REACTOME_SIGNALING_BY_INTERLEUKINS", "REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS",
                   "REACTOME_HSF1_ACTIVATION", "KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY", "WP_CYTOSOLIC_DNASENSING_PATHWAY", "GOBP_ALPHA_BETA_T_CELL_ACTIVATION",
                   "GOBP_T_CELL_DIFFERENTIATION", "GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING", "REACTOME_PI3K_AKT_SIGNALING_IN_CANCER",
                   "PID_CXCR4_PATHWAY", "GOBP_TRNA_MODIFICATION")) %>% 
  select(ID, setSize, NES, p.adjust) %>% 
  mutate(log10p = -log10(p.adjust)) %>% 
  arrange(NES)
GSEA_dot$ID <- factor(GSEA_dot$ID, levels = GSEA_dot$ID, labels = c("HSF1_activation", "Response to heat stress", "tRNA modification",
                                                                    "Protein folding", "T cell differentiation", "CXCR4 pathway",
                                                                    "T cell activation", "Interleukin signaling", "PI3K/AKT signaling",
                                                                    "Cytokine signaling", "Cytosolic DNA sensing", "RIG-I-liike receptor signaling"))

ggplot(GSEA_dot, aes(x = NES, y = ID)) +
  geom_point(aes(size = setSize, fill = log10p), color = "black", shape = 21, stroke = 0.5) +
  geom_vline(xintercept = 0, color = "#1B1919", linetype = "dashed") +
  scale_fill_gradient2(low = "#FCFCFC", mid = "#608AC8", high = "#112641", midpoint = 3) +
  theme_niwot() + 
  theme(panel.border = element_rect(color = "black", size = 2, fill = NA)) +
  labs(x = "NES", y = "")
ggsave(filename = "results/4_funtional_analysis/2_GSEA_dotplot.pdf", height = 4.5, width = 6)

### 3. Signature analysis
rm(list = ls()[-which(ls() %in% c("alldata", "CD8T"))])
gc()
library(GSEABase)
## homologous gene conversion
Homologous_gene <- read_tsv(file = "~/Bioinformatics/Programs/2023/02_yanglichao_scTransSpecies/mart_export.txt") %>% dplyr::select(1:9)

immu_gmt <- getGmt("~/Bioinformatics/Programs/2022/8_public_database/data/DIY_gmt/hsa_Immune_response.gmt")
genesets <- lapply(names(immu_gmt), function(x) {
  gene <- immu_gmt[[x]]@geneIds
  tmp <- Homologous_gene %>% filter(`Gene name` %in% gene) %>% pull(`Mouse gene name`) %>% 
    na.omit() %>% unique()
  tmp <- intersect(tmp, rownames(CD8T))
  return(tmp)
})
names(genesets) <- names(immu_gmt)
genesets <- genesets[c(1,5,10,12:17,19:23,29,39:46,48:54,57:67,75:77,87,93:98)]

library(UCell)
score_c <- AddModuleScore_UCell(CD8T, features = genesets, name = NULL, ncores = 12) %>% 
  .@meta.data %>% dplyr::select(c("group", "celltype1", names(genesets)))
write_csv(score_c, file = "results/4_funtional_analysis/3_signature_Ucell.csv")

library(AUCell)
cells_rankings <- AUCell_buildRankings(GetAssayData(CD8T, layer = "data"))
score_a <- AUCell_calcAUC(geneSets = genesets, rankings = cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05) %>% 
  getAUC() %>% t() %>% as.data.frame() %>% cbind(CD8T@meta.data[,c("group", "celltype1")], .)
write_csv(score_a, file = "results/4_funtional_analysis/3_signature_AUcell.csv")

library(ggsignif)
library(FSA) # Dunn's test
signature_meta <- data.frame(Signatures = names(genesets), Labels = paste0("S", 1:length(genesets)))
data <- score_c %>% dplyr::select(-1) %>% filter(celltype1 %in% c("C1", "C3", "C4"))
colnames(data) <- c("celltype1", signature_meta$Labels)

# data transforming
data <- data %>% 
  pivot_longer(cols = -celltype1, names_to = "Signature", values_to = "Scores")
data$celltype1 <- factor(data$celltype1, levels = c("C1", "C3", "C4"))
data$Signature <- factor(data$Signature, levels = paste0("S", 1:length(genesets)))

# Kruskal-Wallis test
kruskal_results <- data %>%
  group_by(Signature) %>%
  summarise(p.value = kruskal.test(Scores ~ celltype1)$p.value)
# Dunn's 检验
dunn_results <- data %>%
  group_by(Signature) %>%
  do(dunnTest(Scores ~ celltype1, data = ., method = "bonferroni")$res) %>%
  filter(P.adj < 0.05)
dunn_results <- dunn_results %>% 
  filter(Comparison == "C3 - C4") %>% 
  left_join(., signature_meta, by = c("Signature" = "Labels"))

df1 <- data %>% filter(Signature %in% c("S1", "S2", "S12"))
ggplot(df1, aes(x = Signature, y = Scores, fill = celltype1)) +
  geom_violin(trim = F, width = 0.9, position = position_dodge(0.8), alpha = 0.5) +
  geom_boxplot(width = 0.08, position = position_dodge(0.8), outlier.shape = NA) +
  theme_niwot() + theme(legend.position = "top") +
  scale_fill_manual(values = c("C1" = "#E0E0E0", "C3" = "#DC0000", "C4" = "#3C5488")) +
  labs(x = "", y = "Scores")
ggsave(filename = "results/4_funtional_analysis/3_signature1.pdf", width = 6, height = 4)

df2 <- data %>% filter(Signature %in% c("S33", "S35", "S38"))
ggplot(df2, aes(x = Signature, y = Scores, fill = celltype1)) +
  geom_violin(trim = F, width = 0.8, position = position_dodge(0.8), alpha = 0.5) +
  geom_boxplot(width = 0.08, position = position_dodge(0.8), outlier.shape = NA) +
  theme_niwot() + theme(legend.position = "top") +
  scale_fill_manual(values = c("C1" = "#E0E0E0", "C3" = "#DC0000", "C4" = "#3C5488")) +
  labs(x = "", y = "Scores")
ggsave(filename = "results/4_funtional_analysis/3_signature2.pdf", width = 6, height = 4)

### 4. density analysis
rm(list = ls()[-which(ls() %in% c("alldata", "CD8T"))])
gc()
umap_coordinations <- CD8T@reductions[["umap.scanorama"]]@cell.embeddings %>% as.data.frame() %>% cbind(CD8T@meta.data[,c("group", "celltype1")], .)
colnames(umap_coordinations)[3:4] <- c("X", "Y")
write_csv(umap_coordinations, file = "results/4_funtional_analysis/4_umap_scanorama_coordinations.csv")


### Trajectory analysis

CD8T_1 <- subset(CD8T, subset = celltype1 != "C6")
CD8T_1$celltype1 <- factor(CD8T_1$celltype1)
CD8T_1 <- RunSlingshot(srt = CD8T_1, group.by = "celltype1", reduction = "umap.scanorama")
CellDimPlot(CD8T_1, group.by = "celltype1", reduction = "umap.scanorama", lineages = paste0("Lineage", 1:2), lineages_span = 0.1, palcolor = colors[6:1])
ggsave(filename = "results/4_funtional_analysis/5_CD8T_lineage.pdf", height = 5, width = 7)

{
  ## slingshot
  library(slingshot)
  library(SingleCellExperiment)
  CD8T.sce <- as.SingleCellExperiment(CD8T_1)
  sim <- slingshot(CD8T.sce, clusterLabels = "celltype1", reducedDim = Embeddings(CD8T_1, reduction = "scanorama"))
  
  pdf(file = "results/4_funtional_analysis/5_test.pdf")
  par(mar = c(5.1, 4.1, 4.1, 14), xpd = T)
  plot(reducedDims(sim)$UMAP.SCANORAMA, col = colors[6:1][sim$celltype1], pch = 16, asp = 1)
  lines(SlingshotDataSet(sim), lwd = 2, type = "lineages", col = "black")
  dev.off()
}

### 5. Monocle3
library(monocle3)

## 数据准备
data <- GetAssayData(CD8T_1, layer = "counts")
cell_meta <- CD8T_1@meta.data
gene_anno <- data.frame(gene_short_name = rownames(data))
rownames(gene_anno) = rownames(data)

## 创建cds对象
cds <- new_cell_data_set(data, cell_metadata = cell_meta, gene_metadata = gene_anno)

## 标准化和PCA降维
cds <- preprocess_cds(cds, num_dim = 50)

## UMAP降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")

## 导入原有reduction embeddings
cds@int_colData$reducedDims$UMAP.original <- cds@int_colData$reducedDims$UMAP # 备份原有
cds.embed <- cds@int_colData$reducedDims$UMAP
scanorama.embed <- Embeddings(CD8T_1, reduction = "umap.scanorama")
scanorama.embed <- scanorama.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- scanorama.embed
plot_cells(cds, reduction_method = "UMAP", color_cells_by = "celltype1")

## 分区
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

## 构建细胞轨迹

# 轨迹学习
cds <- learn_graph(cds)
plot_cells(cds, color_cells_by = "celltype1", label_groups_by_cluster = T, label_leaves = T, label_branch_points = T, graph_label_size = 1.5)

# 选择root
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = F, label_roots = T, label_leaves = F, label_branch_points = T)

## 将pseudotime信息添加至原数据
identical(rownames(cell_meta), names(cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]))
CD8T_1$pseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
FeatureDimPlot(CD8T_1, features = "pseudotime", reduction = "umap.scanorama", palette = "YlOrRd")
ggsave(filename = "results/4_funtional_analysis/5_CD8T_pseudotime.pdf", height = 5, width = 5)

## 手工绘制root, branch, leaves
source("~/Bioinformatics/Programs/2022/8_public_database/codes/ggplot_theme.R", echo=TRUE)
data_df <- cbind(Embeddings(CD8T_1, reduction = "umap.scanorama"), select(CD8T_1@meta.data, celltype1, pseudotime))
colnames(data_df)[1:2] <- c("UMAP_1", "UMAP_2")

## data processing
{
  # Define the branch_nodes function
  branch_nodes <- function(cds, reduction_method) {
    # Extract the principal graph
    dp_mst <- cds@principal_graph[[reduction_method]]
    
    # Identify branch nodes (nodes with more than two connections)
    branch_nodes <- which(igraph::degree(dp_mst) > 2)
    
    return(branch_nodes)
  }
  
  # Define the leaf_nodes function
  leaf_nodes <- function(cds, reduction_method) {
    # Extract the principal graph
    dp_mst <- cds@principal_graph[[reduction_method]]
    
    # Identify leaf nodes (nodes with only one connection)
    leaf_nodes <- which(igraph::degree(dp_mst) == 1)
    
    return(leaf_nodes)
  }
  
  # Define the root_nodes function
  root_nodes <- function(cds, reduction_method) {
    # Extract the principal graph
    dp_mst <- cds@principal_graph[[reduction_method]]
    
    # Identify root nodes (nodes with zero or one connection)
    root_nodes <- which(igraph::degree(dp_mst) <= 1)
    
    return(root_nodes)
  }
  
  # tables
  ica_space_df <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst) %>% as.data.frame() %>% 
    dplyr::select(prin_graph_dim_1 = umapscanorama_1, prin_graph_dim_2 = umapscanorama_2) %>% 
    dplyr::mutate(sample_name = rownames(.), sample_state = rownames(.))
  dp_mst <- cds@principal_graph[["UMAP"]]
  edge_df <- dp_mst %>% igraph::as_data_frame() %>% 
    dplyr::select(source = "from", target = "to") %>% 
    dplyr::left_join(ica_space_df %>% dplyr::select(source = "sample_name", 
                                                    source_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                    source_prin_graph_dim_2 = "prin_graph_dim_2"), 
                     by = "source") %>% 
    dplyr::left_join(ica_space_df %>% dplyr::select(target = "sample_name", 
                                                    target_prin_graph_dim_1 = "prin_graph_dim_1", 
                                                    target_prin_graph_dim_2 = "prin_graph_dim_2"), 
                     by = "target")
}

## plot
p <- ggplot(data_df, aes(x = UMAP_1, y = UMAP_2, color = pseudotime)) +
  geom_point(size = 1, alpha = 1) +
  scale_color_continuous_sequential(palette = "YlGnBu") + 
  theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white")) +
  theme(
    legend.title = element_text(size = 12, face = "bold", vjust = 1),
    legend.text = element_text(size = 10),
    legend.direction = "horizontal",
    legend.position = c(0.7,0.9),
  ) +
  coord_fixed() +
  labs(color = "Pseudotime")

# ## celltype1 legend
# theme(
#   legend.title = element_blank(), #去掉legend.title 
#   legend.key=element_rect(fill='white'), #
#   legend.text = element_text(size=20), #设置legend标签的大小
#   legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
#   guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小

## XY轴annotation
p <- p + 
  geom_segment(aes(x = min(data_df$UMAP_1) , y = min(data_df$UMAP_2),
                   xend = min(data_df$UMAP_1) + 1.5, yend = min(data_df$UMAP_2)),
               colour = "black", size = 1, arrow = arrow(length = unit(0.3,"cm"))) + 
  geom_segment(aes(x = min(data_df$UMAP_1)  , y = min(data_df$UMAP_2),
                   xend = min(data_df$UMAP_1), yend = min(data_df$UMAP_2) + 1.5),
               colour = "black", size = 1, arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(data_df$UMAP_1) + 0.5, y = min(data_df$UMAP_2) - 0.3, label = "UMAP_1",
           color = "black", size = 4, fontface = "bold") + 
  annotate("text", x = min(data_df$UMAP_1) - 0.3, y = min(data_df$UMAP_2) + 0.5, label = "UMAP_2",
           color = "black", size = 4, fontface = "bold", angle = 90)

## show_trajectory_graph
trajectory_graph_segment_size <- 0.75
graph_label_size <- 3
p <- p + geom_segment(data = edge_df, 
                 aes_string(x = "source_prin_graph_dim_1", 
                            y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                            yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size, 
                 color = I("grey28"), linetype = "solid", 
                 na.rm = TRUE)

## label_roots
mst_root_nodes <- root_nodes(cds, reduction_method = "UMAP")
# set root nodes: index 3
root_df <- ica_space_df %>% 
  dplyr::slice(match(names(mst_root_nodes), sample_name)) %>% 
  dplyr::mutate(root_idx = seq_len(dplyr::n()))
root_df <- root_df %>% filter(sample_name == "Y_50")
root_df$root_reidx <- "1"
p <- p + geom_point(data = root_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                    color = "black", fill = "white", size = I(graph_label_size * 1.5), na.rm = TRUE) + 
  geom_text(data = root_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "root_reidx"), 
            size = I(graph_label_size), color = "black", na.rm = TRUE)

## label_leaves
mst_leaf_nodes <- leaf_nodes(cds, reduction_method = "UMAP")
# set leave nodes
leaf_df <- ica_space_df %>% 
  dplyr::slice(match(names(mst_leaf_nodes), sample_name)) %>% 
  dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
leaf_df <- leaf_df %>% filter(sample_name != "Y_50")
leaf_df$leaf_reidx <- c("1", "2")
p <- p + geom_point(data = leaf_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
                    color = "black", fill = "lightgray", size = I(graph_label_size * 1.5), na.rm = TRUE) + 
  geom_text(data = leaf_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "leaf_reidx"), 
            size = I(graph_label_size), color = "black", na.rm = TRUE)

## label_branch_point
mst_branch_nodes <- branch_nodes(cds, reduction_method = "UMAP")
# set branch_nodes
branch_point_df <- ica_space_df %>% dplyr::slice(match(names(mst_branch_nodes), sample_name)) %>% 
  dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
branch_point_df$branch_point_reidx <- c("4", "7", "1", "3", "2", "5", "6", "8")
p <- p + geom_point(data = branch_point_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), 
                    shape = 21, stroke = I(trajectory_graph_segment_size),  color = "white", fill = "black", 
                    size = I(graph_label_size * 1.5), na.rm = TRUE) + 
  geom_text(data = branch_point_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "branch_point_reidx"), 
            size = I(graph_label_size), color = "white", na.rm = TRUE)
print(p)

FeatureDimPlot(CD8T_1, features = "pseudotime", reduction = "umap.scanorama", palette = "YlOrRd") +
  geom_segment(data = edge_df, 
               aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), size = trajectory_graph_segment_size, 
               color = I("grey28"), linetype = "solid", 
               na.rm = TRUE) +
  geom_point(data = root_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
             color = "black", fill = "white", size = I(graph_label_size * 1.2), na.rm = TRUE) + 
  geom_text(data = root_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "root_reidx"), 
            size = I(graph_label_size), color = "black", na.rm = TRUE) +
  geom_point(data = leaf_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), shape = 21, stroke = I(trajectory_graph_segment_size), 
             color = "black", fill = "lightgray", size = I(graph_label_size * 1.2), na.rm = TRUE) + 
  geom_text(data = leaf_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "leaf_reidx"), 
            size = I(graph_label_size), color = "black", na.rm = TRUE) +
  geom_point(data = branch_point_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2"), 
             shape = 21, stroke = I(trajectory_graph_segment_size),  color = "white", fill = "black", 
             size = I(graph_label_size * 1.2), na.rm = TRUE) + 
  geom_text(data = branch_point_df, aes_string(x = "prin_graph_dim_1", y = "prin_graph_dim_2", label = "branch_point_reidx"), 
            size = I(graph_label_size), color = "white", na.rm = TRUE)
ggsave(filename = "results/4_funtional_analysis/5_CD8T_pseudotime2.pdf", height = 5, width = 5)

### 6. SCENIC analysis
## convert seuratobj to anndata
adata <- srt_to_adata_v5(srt = alldata, 
                         meta_to_use = c("orig.ident", "nCount_RNA", "nFeature_RNA", "sample", "group", "main_anno", "fine_anno"),
                         reduction_to_use = c("scanorama", "umap.scanorama"))
## Save data
adata$write_h5ad("results/5_SCENIC_analysis/alldata_to_anndata_raw.h5ad")

### 7. TCR analysis

## clonal expansion
Tcell <- rbind(CD8T@meta.data %>% rownames_to_column("CellID") %>% 
                 mutate(celltype2 = paste0("CD8T_", celltype1)) %>% select(CellID, fine_anno = celltype2),
               subset(alldata, subset = main_anno == "CD4T") %>% .@meta.data %>% 
                 rownames_to_column("CellID") %>% select(CellID, fine_anno))
write_csv(Tcell, file = "results/6_TCR_analysis/Tcell_to_analysis.csv")

TCR_profile <- read_csv(file = "results/6_TCR_analysis/TCR_profile.csv") %>% 
  mutate(CellID = str_sub(CellID, 1, -3))

CD8T_TCR <- CD8T[, TCR_profile$CellID]
CD8T_TCR$clonal_expansion <- factor(TCR_profile %>% column_to_rownames("CellID") %>% .[colnames(CD8T_TCR),"airr:clonal_expansion"])
CellDimPlot(CD8T_TCR, group.by = "clonal_expansion", reduction = "umap.scanorama", split.by = "group", palette = "inferno")
ggsave(filename = "results/6_TCR_analysis/TCR_clonal_expansion_umap.pdf", height = 5, width = 10)

data <- CD8T_TCR@meta.data %>% select(group, celltype1, clonal_expansion)
data$clonal_expansion <- factor(data$clonal_expansion, levels = c("> 1000", "<= 1000", "<= 100", "<= 10", "<= 1"))
# 将clonal_expansion转换为数值型
data$clonal_expansion2 <- as.numeric(data$clonal_expansion)
# 转换clonal_expansion为百分比
data <- transform(data, percent = clonal_expansion2 / ave(clonal_expansion2, celltype1, FUN = sum) * 100)
# 绘制堆叠百分比图
p <- ggplot(data, aes(x = celltype1, y = percent, fill = clonal_expansion)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  facet_grid(. ~ group) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Cell Type", y = "Percentage of Clonal Expansion", fill = "Clonal Expansion") +
  scale_fill_manual(values = c("#FCFFA4", "#F98C09", "#BB3654", "#560F6D", "#000004")) +
  theme_niwot() +
  theme(panel.background = element_rect(color = "black"))
save_plot(filename = "results/6_TCR_analysis/TCR_clonal_expasion_bar_stack.pdf", plot = p, height = 5, width = 6, dpi = 300)

## shared clone
weightedMtx <- read_csv(file = "results/6_TCR_analysis/TCR_weighted_matrix.csv")
weightedMtx <- column_to_rownames(weightedMtx, "gex:fine_anno_by_group")
sum(weightedMtx)
shared_clone <- c()
for (j in colnames(weightedMtx)) {
  if (sum(weightedMtx[,j] > 0) > 1) {
    shared_clone <- c(shared_clone, j)
  }
}
shared_stat <- data.frame(shared = sapply(rownames(weightedMtx), function(i) sum(weightedMtx[i, shared_clone])),
                          total = sapply(rownames(weightedMtx), function(i) sum(weightedMtx[i, ])),
                          percent = sapply(rownames(weightedMtx), function(i) sum(weightedMtx[i, shared_clone]) / sum(weightedMtx[i,]) * 100),
                          row.names = rownames(weightedMtx)) %>% 
  rownames_to_column("fine_anno_by_group") %>% 
  separate(col = "fine_anno_by_group", into = c("group", "split", "fine_anno"), sep = c(3,4)) %>% 
  select(-split)
write_csv(shared_stat, file = "results/6_TCR_analysis/TCR_shared_stat.csv")
shared_stat <- shared_stat %>%  filter(grepl("^CD8T", fine_anno)) %>% 
  rename(yes = percent) %>% mutate(no = 100 - yes) %>% 
  pivot_longer(cols = c("yes", "no"), names_to = "shared", values_to = "percentage") %>% 
  mutate(label = ifelse(shared == "yes", paste0(round(percentage, 1), "%"), ""))
shared_stat$group <- factor(shared_stat$group, levels = c("ctr", "cko"))
shared_stat$shared <- factor(shared_stat$shared, levels = c("yes", "no"))

ggplot(data = shared_stat, aes(x = "", y = percentage, fill = shared)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  geom_text(aes(label = label), position = position_stack(0.5), size = 3) +
  coord_polar(theta = "y", start = 0, direction = -1) + 
  facet_grid(group ~ fine_anno) +
  theme_void() +
  scale_fill_manual(values = c("#A30034", "#FCF5F2"))
ggsave(filename = "results/6_TCR_analysis/TCR_shared_clone.pdf", width = 8, height = 3)

## Chi-square
C2 <- matrix(c(138, 68, 230, 283), nrow=2, ncol=2)
chisq.test(C2)

C3 <- matrix(c(1310, 185, 97, 26), nrow=2, ncol=2)
chisq.test(C3)

## distance
distance <- read_csv(file = "results/6_TCR_analysis/TCR_distance.csv")
df <- distance %>% column_to_rownames("gex:fine_anno_by_group") %>% as.matrix()
group <- c("ctr", "cko")
celltype <- paste0("CD8T_C", 1:6)
meta <- expand_grid(group, celltype) %>% 
  mutate(cellname = paste(group, celltype, sep = "_"))
df <- df[meta$cellname, meta$cellname]

col_fun <- colorRamp2(c(0, 0.2, 0.4), c("#440459", "#1F928B", "#F8E620"))
col_anno <- HeatmapAnnotation(Group = meta$group, `Cell type` = meta$celltype,
                              col = list(
                                Group = c("ctr" = "#B0A875", "cko" = "#AF7EC0"),
                                `Cell type` = c("CD8T_C1" = "#E31A1C", "CD8T_C2" = "#FB9A99", "CD8T_C3" = "#33A02C",
                                                "CD8T_C4" = "#B2DF8A", "CD8T_C5" = "#1F78B4", "CD8T_C6" = "#A6CEE3")
                              )
)
row_anno <- rowAnnotation(Group = meta$group, `Cell type` = meta$celltype,
                          col = list(
                            Group = c("ctr" = "#B0A875", "cko" = "#AF7EC0"),
                            `Cell type` = c("CD8T_C1" = "#E31A1C", "CD8T_C2" = "#FB9A99", "CD8T_C3" = "#33A02C",
                                            "CD8T_C4" = "#B2DF8A", "CD8T_C5" = "#1F78B4", "CD8T_C6" = "#A6CEE3")
                          ),
                          show_legend = c(FALSE, FALSE),
                          show_annotation_name = F)
pdf(file = "results/6_TCR_analysis/TCR_distance_overlap.pdf", width = 6, height = 5)
Heatmap(df, na_col = "grey95", col = col_fun, 
        top_annotation = col_anno, left_annotation = row_anno,
        show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F,
        width = ncol(df)*unit(10, "mm"), height = nrow(df)*unit(10, "mm"),
        rect_gp = gpar(col = "white", lwd = 2),
        heatmap_legend_param = list(title = "1-Distance"))
dev.off()

## relative distance
df1 <- df[meta$cellname[1:6], meta$cellname[1:6]]
df2 <- df[meta$cellname[7:12], meta$cellname[7:12]]

df_r <- df2 - df1
colnames(df_r) <- celltype
rownames(df_r) <- celltype

col_fun2 <- colorRamp2(c(-0.07, 0, 0.05), c("#023FA5", "white", "#A30034"))
col_anno2 <- HeatmapAnnotation(`Cell type` = celltype,
                              col = list(
                                `Cell type` = c("CD8T_C1" = "#E31A1C", "CD8T_C2" = "#FB9A99", "CD8T_C3" = "#33A02C",
                                                "CD8T_C4" = "#B2DF8A", "CD8T_C5" = "#1F78B4", "CD8T_C6" = "#A6CEE3")
                              ),
                              show_annotation_name = F
)
row_anno2 <- rowAnnotation(`Cell type` = celltype,
                          col = list(
                            `Cell type` = c("CD8T_C1" = "#E31A1C", "CD8T_C2" = "#FB9A99", "CD8T_C3" = "#33A02C",
                                            "CD8T_C4" = "#B2DF8A", "CD8T_C5" = "#1F78B4", "CD8T_C6" = "#A6CEE3")
                          ),
                          show_legend = c(FALSE, FALSE),
                          show_annotation_name = F)
pdf(file = "results/6_TCR_analysis/TCR_relative_distance.pdf", width = 4, height = 3)
Heatmap(df_r, na_col = "grey95", col = col_fun2, 
        top_annotation = col_anno2, left_annotation = row_anno2,
        show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F,
        width = ncol(df_r)*unit(5, "mm"), height = nrow(df_r)*unit(5, "mm"),
        rect_gp = gpar(col = "black", lwd = 1),
        heatmap_legend_param = list(title = "Relative 1-Distance"))
dev.off()


#### CD4T ----
rm(list = ls()[-which(ls() %in% c("alldata"))])
gc()

CD4T <- subset(alldata, subset = main_anno == "CD4T")
CellDimPlot(CD4T, group.by = "fine_anno", reduction = "umap.scanorama", label = T)
selected_cells <- Embeddings(CD4T, reduction = "umap.scanorama") %>% as.data.frame() %>% 
  filter(umapscanorama_2 <= -2 & umapscanorama_1 <= 1) %>% rownames()
CD4T <- subset(CD4T, cells = selected_cells)
CellDimPlot(CD4T, group.by = "fine_anno", reduction = "umap.scanorama", label = T)
CD4T$celltype1 <- factor(CD4T$fine_anno)
levels(CD4T$celltype1) <- paste0("C", 1:6)
Idents(CD4T) <- "celltype1"
CD4_markers <- FindAllMarkers(CD4T)
write_csv(CD4_markers, file = "results/4_funtional_analysis/6_CD4_markers.csv")
saveRDS(CD4T, file = "results/8_rebuttal/6_CD4T.rds")

### 6. UMAP
colors <- brewer.pal(12, "Paired")
CellDimPlot(CD4T, group.by = "celltype1", reduction = "umap.scanorama", label = T, palcolor = colors[6:1])
ggsave(filename = "results/4_funtional_analysis/6_CD4T_cluster.pdf", height = 5, width = 5)
CellDimPlot(CD4T, group.by = "celltype1", reduction = "umap.scanorama", split.by = "group", label = T, palcolor = colors[6:1])
ggsave(filename = "results/4_funtional_analysis/6_CD4T_cluster_by_group.pdf", height = 5, width = 10)

### 6. cell proportion
CellStatPlot(CD4T, stat.by = "celltype1", group.by = "group", plot_type = "trend", label = T, palcolor = colors[6:1])
ggsave(filename = "results/4_funtional_analysis/6_CD4T_proportion.pdf", width = 3.5, height = 5)

### 6. marker gene umap
VlnPlot(CD4T, features = c("Foxp3", "Il2ra", "Sell", "Cd44", "Tbx21", "Gata3", "Ifng", "Cd40lg", "Tnfsf9"))
FeatureDimPlot_v5(CD4T, features = c("Foxp3", "Il2ra", "Sell", "Cd44", "Tbx21", "Gata3", "Ifng", "Cd40lg", "Tnfsf9"), reduction = "umap.scanorama")
ggsave(filename = "results/4_funtional_analysis/6_CD4T_gene_umap.pdf", height = 8, width = 10)

FeatureDimPlot_v5(CD4T, features = c("Hspa1a", "Hspa1b", "Hsph1", "Dnajb1", "Klf2", "Foxp1"), reduction = "umap.scanorama", split.by = "group")
ggsave(filename = "results/4_funtional_analysis/6_CD4T_gene_umap_hsp.pdf", height = 8, width = 12)

FeatureStatPlot_v5(subset(CD4T, subset = celltype1 %in% c("C3", "C4", "C5")), stat.by = c("Ctcf", "Foxp3", "Il2ra", "Tbx21", "Ifng", "Cd40lg"), 
                   group.by = "celltype1", split.by = "group", comparisons = TRUE)
expr <- AverageExpression(subset(CD4T, subset = celltype1 %in% c("C3")), assays = "RNA", 
                          features = c("Foxp3", "Il2ra", "Cst7", "Nkg7", "Gzma", "Prf1", "Tbx21", "Ccl4", "Ccl3", "Ccl5"), group.by = "group")$RNA %>% as.matrix()
ggsave(filename = "results/4_funtional_analysis/6_CD4T_gene_violinplot.pdf", height = 6, width = 10)

### 7. Dotplot
genes2plot <- c("Cd44", "Sell", "Il7r", "Foxp3", "Il2ra", "Tnfrsf18", "Tbx21", "Gata3", "Rorc", "Ifng", "Lta", "Il10", "Il17a", "Cd40lg", "Stmn1")
DotExpr <- AverageExpression(CD4T, assays = "RNA", features = genes2plot, group.by = "celltype1")$RNA %>% as.matrix()
DotExpr <- log2(DotExpr[genes2plot,] + 1)
Idents(CD4T) <- "celltype1"
markers1 <- FindAllMarkers(CD4T, features = genes2plot, logfc.threshold = 0, min.pct = 0)
DotData <- as.data.frame(DotExpr) %>% rownames_to_column("gene") %>% 
  gather(key = "cluster", value = "expression", -gene) %>% 
  left_join(., select(markers1, c(gene, cluster, pct.1)), by = c("gene", "cluster"))
DotData$pct.1 <- ifelse(is.na(DotData$pct.1), 0, DotData$pct.1)
DotData$gene <- factor(DotData$gene, levels = rev(genes2plot))

ggplot(DotData, aes(x = cluster, y = gene)) +
  geom_point(aes(size = pct.1, fill = expression), color = "black", shape = 21, stroke = 0.5) +
  scale_fill_gradient2(low = "#F6F6F6", mid = "#FF888E", high = "#A01021", midpoint = 2.5) +
  theme_niwot() +
  labs(x = "", y = "")
ggsave(filename = "results/4_funtional_analysis/7_CD4T_dotplot.pdf", height = 4.5, width = 4.5)

### 8. Differential expression analysis
cko <- subset(CD4T, subset = group == "cko")
Idents(cko) <- "celltype1"
cko_C4vsC5 <- FindMarkers(cko, ident.1 = "C4", ident.2 = "C5", logfc.threshold = 0) %>% rownames_to_column("gene")
write_csv(cko_C4vsC5, file = "results/4_funtional_analysis/8_CD4T_cko_C4vsC5.csv")

C4 <- subset(CD4T, subset = celltype1 == "C4")
Idents(C4) <- "group"
C4_cKOvsctr <- FindMarkers(C4, ident.1 = "cko", ident.2 = "ctr", logfc.threshold = 0) %>% rownames_to_column("gene")
write_csv(C4_cKOvsctr, file = "results/4_funtional_analysis/8_CD4T_C4_cKOvsctr.csv")

### 9. Signature analysis
rm(list = ls()[-which(ls() %in% c("alldata", "CD4T"))])
gc()
library(GSEABase)
## homologous gene conversion
Homologous_gene <- read_tsv(file = "~/Bioinformatics/Programs/2023/02_yanglichao_scTransSpecies/mart_export.txt") %>% dplyr::select(1:9)

immu_gmt <- getGmt("/home/data/qiuli/Bioinformatics/Programs/2022/8_public_database/data/DIY_gmt/hsa_Immune_response.gmt")
genesets <- lapply(names(immu_gmt), function(x) {
  gene <- immu_gmt[[x]]@geneIds
  tmp <- Homologous_gene %>% filter(`Gene name` %in% gene) %>% pull(`Mouse gene name`) %>% 
    na.omit() %>% unique()
  tmp <- intersect(tmp, rownames(CD4T))
  return(tmp)
})
names(genesets) <- names(immu_gmt)
genesets <- genesets[c(3,17,24:25,69:71,92)]

library(UCell)
score_c <- AddModuleScore_UCell(CD8T, features = genesets, name = NULL, ncores = 12) %>% 
  .@meta.data %>% dplyr::select(c("group", "celltype1", names(genesets)))
write_csv(score_c, file = "results/4_funtional_analysis/3_signature_Ucell.csv")

library(AUCell)
cells_rankings <- AUCell_buildRankings(GetAssayData(CD8T, layer = "data"))
score_a <- AUCell_calcAUC(geneSets = genesets, rankings = cells_rankings, aucMaxRank = nrow(cells_rankings) * 0.05) %>% 
  getAUC() %>% t() %>% as.data.frame() %>% cbind(CD8T@meta.data[,c("group", "celltype1")], .)
write_csv(score_a, file = "results/4_funtional_analysis/3_signature_AUcell.csv")

library(ggsignif)
library(FSA) # Dunn's test
signature_meta <- data.frame(Signatures = names(genesets), Labels = paste0("S", 1:length(genesets)))
data <- score_c %>% dplyr::select(-1) %>% filter(celltype1 %in% c("C1", "C3", "C4"))
colnames(data) <- c("celltype1", signature_meta$Labels)

# data transforming
data <- data %>% 
  pivot_longer(cols = -celltype1, names_to = "Signature", values_to = "Scores")
data$celltype1 <- factor(data$celltype1, levels = c("C1", "C3", "C4"))
data$Signature <- factor(data$Signature, levels = paste0("S", 1:length(genesets)))

# Kruskal-Wallis test
kruskal_results <- data %>%
  group_by(Signature) %>%
  summarise(p.value = kruskal.test(Scores ~ celltype1)$p.value)
# Dunn's 检验
dunn_results <- data %>%
  group_by(Signature) %>%
  do(dunnTest(Scores ~ celltype1, data = ., method = "bonferroni")$res) %>%
  filter(P.adj < 0.05)
dunn_results <- dunn_results %>% 
  filter(Comparison == "C3 - C4") %>% 
  left_join(., signature_meta, by = c("Signature" = "Labels"))

df1 <- data %>% filter(Signature %in% c("S1", "S2", "S12"))
ggplot(df1, aes(x = Signature, y = Scores, fill = celltype1)) +
  geom_violin(trim = F, width = 0.9, position = position_dodge(0.8), alpha = 0.5) +
  geom_boxplot(width = 0.08, position = position_dodge(0.8), outlier.shape = NA) +
  theme_niwot() + theme(legend.position = "top") +
  scale_fill_manual(values = c("C1" = "#E0E0E0", "C3" = "#DC0000", "C4" = "#3C5488")) +
  labs(x = "", y = "Scores")
ggsave(filename = "results/4_funtional_analysis/3_signature1.pdf", width = 6, height = 4)

### 10. TCR analysis
rm(list = ls()[-which(ls() %in% c("alldata", "CD4T"))])
gc()

## clonal expansion

TCR_profile <- read_csv(file = "results/6_TCR_analysis/TCR_profile.csv") %>% 
  mutate(CellID = str_sub(CellID, 1, -3))
CD4_TCR_profile <- TCR_profile %>% filter(`gex:main_anno` == "CD4T") %>% 
  mutate(clonal_expansion = cut(`airr:clone_id_size`, 
                                breaks = c(-Inf, 1, 2, 10, 50, Inf),
                                labels = c("<= 1", "<= 2", "<= 10", "<= 50", "> 50"),
                                include.lowest = TRUE))

CD4_TCR <- CD4T[, CD4_TCR_profile$CellID]
CD4_TCR_profile <- CD4_TCR_profile %>% filter(CellID %in% colnames(CD4_TCR))
CD4_TCR$clonal_expansion <- factor(CD4_TCR_profile %>% column_to_rownames("CellID") %>% .[colnames(CD4_TCR),"clonal_expansion"])
CellDimPlot(CD4_TCR, group.by = "clonal_expansion", reduction = "umap.scanorama", split.by = "group", palette = "inferno")
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_clonal_expansion_umap.pdf", height = 5, width = 10)

data <- CD4_TCR@meta.data %>% dplyr::select(group, celltype1, clonal_expansion)
data$clonal_expansion <- factor(data$clonal_expansion, levels = c("> 50", "<= 50", "<= 10", "<= 2", "<= 1"))
# 将clonal_expansion转换为数值型
data$clonal_expansion2 <- as.numeric(data$clonal_expansion)
# 转换clonal_expansion为百分比
data <- transform(data, percent = clonal_expansion2 / ave(clonal_expansion2, celltype1, FUN = sum) * 100)
# 绘制堆叠百分比图
p <- ggplot(data, aes(x = celltype1, y = percent, fill = clonal_expansion)) +
  geom_bar(stat = "identity", position = "fill", width = 0.7) +
  facet_grid(. ~ group) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_y_continuous(expand = c(0,0)) +
  labs(x = "Cell Type", y = "Percentage of Clonal Expansion", fill = "Clonal Expansion") +
  scale_fill_manual(values = c("#FCFFA4", "#F98C09", "#BB3654", "#560F6D", "#000004")) +
  theme_niwot() +
  theme(panel.background = element_rect(color = "black"))
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_clonal_expasion_bar_stack.pdf", plot = p, height = 5, width = 6, dpi = 300)

## shared clone
weightedMtx <- read_csv(file = "results/6_TCR_analysis/TCR_weighted_matrix.csv")
weightedMtx <- column_to_rownames(weightedMtx, "gex:fine_anno_by_group")
sum(weightedMtx)
shared_clone <- c()
for (j in colnames(weightedMtx)) {
  if (sum(weightedMtx[,j] > 0) > 1) {
    shared_clone <- c(shared_clone, j)
  }
}
shared_stat <- sapply(rownames(weightedMtx), function(i) {
  percent <- sum(weightedMtx[i, shared_clone]) / sum(weightedMtx[i,]) * 100
}) %>% 
  as.data.frame() %>% dplyr::rename(percent = ".") %>% rownames_to_column("fine_anno_by_group") %>% 
  separate(col = "fine_anno_by_group", into = c("group", "split", "fine_anno"), sep = c(3,4)) %>% 
  dplyr::select(-split)
write_csv(shared_stat, file = "results/6_TCR_analysis/TCR_shared_stat.csv")
shared_stat <- shared_stat %>%  filter(grepl("^CD4T", fine_anno)) %>% 
  dplyr::rename(yes = percent) %>% mutate(no = 100 - yes) %>% 
  pivot_longer(cols = c("yes", "no"), names_to = "shared", values_to = "percentage") %>% 
  mutate(label = ifelse(shared == "yes", paste0(round(percentage, 1), "%"), ""))
shared_stat$group <- factor(shared_stat$group, levels = c("ctr", "cko"))
shared_stat$shared <- factor(shared_stat$shared, levels = c("yes", "no"))
shared_stat <- shared_stat %>% 
  mutate(celltype1 = factor(shared_stat$fine_anno, levels = c(paste0("CD4T_C", 0:5)), labels = c(paste0("C", 1:6))))

ggplot(data = shared_stat, aes(x = "", y = percentage, fill = shared)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  geom_text(aes(label = label), position = position_stack(0.5), size = 3) +
  coord_polar(theta = "y", start = 0, direction = -1) + 
  facet_grid(group ~ celltype1) +
  theme_void() +
  scale_fill_manual(values = c("#A30034", "#FCF5F2"))
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_shared_clone.pdf", width = 8, height = 3)

## distance
distance <- read_csv(file = "results/6_TCR_analysis/CD4_TCR_distance.csv")
df <- distance %>% column_to_rownames("gex:fine_anno_by_group") %>% as.matrix()
group <- c("ctr", "cko")
celltype <- paste0("CD4T_C", 0:5)
meta <- expand_grid(group, celltype) %>% 
  mutate(cellname = paste(group, celltype, sep = "_"))
df <- df[meta$cellname, meta$cellname]

col_fun <- colorRamp2(c(0, 0.12, 0.24), c("#440459", "#1F928B", "#F8E620"))
col_anno <- HeatmapAnnotation(Group = meta$group, `Cell type` = meta$celltype,
                              col = list(
                                Group = c("ctr" = "#B0A875", "cko" = "#AF7EC0"),
                                `Cell type` = c("CD4T_C0" = "#E31A1C", "CD4T_C1" = "#FB9A99", "CD4T_C2" = "#33A02C",
                                                "CD4T_C3" = "#B2DF8A", "CD4T_C4" = "#1F78B4", "CD4T_C5" = "#A6CEE3")
                              )
)
row_anno <- rowAnnotation(Group = meta$group, `Cell type` = meta$celltype,
                          col = list(
                            Group = c("ctr" = "#B0A875", "cko" = "#AF7EC0"),
                            `Cell type` = c("CD4T_C0" = "#E31A1C", "CD4T_C1" = "#FB9A99", "CD4T_C2" = "#33A02C",
                                            "CD4T_C3" = "#B2DF8A", "CD4T_C4" = "#1F78B4", "CD4T_C5" = "#A6CEE3")
                          ),
                          show_legend = c(FALSE, FALSE),
                          show_annotation_name = F)
pdf(file = "results/6_TCR_analysis/CD4_TCR_distance_overlap.pdf", width = 6, height = 5)
Heatmap(df, na_col = "grey95", col = col_fun, 
        top_annotation = col_anno, left_annotation = row_anno,
        show_row_names = F, show_column_names = F, cluster_rows = F, cluster_columns = F,
        width = ncol(df)*unit(10, "mm"), height = nrow(df)*unit(10, "mm"),
        rect_gp = gpar(col = "white", lwd = 2),
        heatmap_legend_param = list(title = "1-Distance"))
dev.off()

# clone的mode转code函数
mode2code <- function(mode) {
  code <- numeric(6)  # 创建长度为6的数值向量，初始值为0
  
  for (i in 1:6) {
    if (grepl(paste0("C", i), mode)) {
      code[i] <- 1
    } else {
      code[i] <- 0
    }
  }
  
  paste(code, collapse = "")  # 将数值向量转换为字符串并返回
}

{
  ## component of CD4T_ctr
  
  CD4_TCR_ctr_profile <- CD4_TCR_profile %>% filter(`gex:group` == "ctr") %>% 
    dplyr::select(CellID, celltype1 = `gex:fine_anno`, clone_id = `airr:clone_id`, clone_id_size = `airr:clone_id_size`)
  CD4_TCR_ctr_profile$celltype1 <- factor(CD4_TCR_ctr_profile$celltype1, levels = paste0("CD4T_C", 0:5), labels = paste0("C", 1:6))
  
  CD4_TCR_ctr_profile <- lapply(unique(CD4_TCR_ctr_profile$clone_id), function(i) {
    clones <- CD4_TCR_ctr_profile %>% filter(clone_id == i)
    mode <- clones$celltype1 %>% unique() %>% sort() %>% paste(collapse = "-")
    clones$shared_mode <- mode
    clones$code <- mode2code(mode)
    return(clones)
  }) %>% bind_rows()
  
  CD4_TCR_ctr_list <- lapply(paste0("C", 1:6), function(i) {
    CD4_TCR_ctr_profile %>% filter(celltype1 == i) %>% pull(clone_id)
  })
  names(CD4_TCR_ctr_list) <- paste0("C", 1:6)
  
  m_ctr = make_comb_mat(CD4_TCR_ctr_list, top_n_sets = 6)
  ss <- set_size(m_ctr)
  cs <- comb_size(m_ctr)
  ss_df <- as.data.frame(table(CD4_TCR_ctr_profile$celltype1)) %>% column_to_rownames(colnames(.)[1])
  ss_vec <- setNames(as.numeric(ss_df[names(ss),1]), names(ss))
  cs_df <- as.data.frame(table(CD4_TCR_ctr_profile$code)) %>% column_to_rownames(colnames(.)[1])
  cs_vec <- setNames(as.numeric(cs_df[names(cs),1]), names(cs))
  ss_mtx <- matrix(ncol = 2, c(ss_vec,ss)) * -1
  rownames(ss_mtx) <- names(ss)
  cs_mtx <- matrix(ncol = 2, c(cs_vec,cs))
  rownames(cs_mtx) <- names(cs)
  ht <- UpSet(m_ctr,
              set_order = 1:6,
              comb_order = order(comb_degree(m_ctr), -cs),
              top_annotation = HeatmapAnnotation(
                "Shared clones" = anno_barplot(cs_mtx,
                                               ylim = c(0, max(cs) * 1.1),
                                               border = F,
                                               gp = gpar(fill = c("black", "white"), col = c("black", "black")),
                                               height = unit(4, "cm"),
                                               beside = T, attach = T),
                annotation_name_side = "left",
                annotation_name_rot = 90
              ),
              left_annotation = rowAnnotation(
                "Clones per cell type" = anno_barplot(ss_mtx,
                                                      baseline = 0,
                                                      axis_param = list(
                                                        at = c(0, -500, -1000),
                                                        labels = c(0, 500, 1000),
                                                        labels_rot = 0),
                                                      border = F,
                                                      gp = gpar(fill = c("black", "white"), col = c("black", "black")),
                                                      width = unit(4, "cm"),
                                                      beside = T, attach = T),
                set_name = anno_text(set_name(m_ctr),
                                     location = 0.5,
                                     just = "center",
                                     width = max_text_width(set_name(m_ctr)) + unit(4, "mm"))
              ),
              right_annotation = NULL,
              show_row_names = F,
              row_names_side = "right",
              width = length(cs)*unit(5, "mm"), height = length(ss)*unit(5, "mm"))
  pdf(file = "results/6_TCR_analysis/CD4_TCR_ctr_clone_profile.pdf", width = 15, height = 5)
  ht <- draw(ht)
  co <- column_order(ht)
  ro <- row_order(ht)
  decorate_annotation("Shared clones", {
    grid.text(cs_mtx[co,1], x = seq_along(cs) - 1/6, y = unit(cs_mtx[co,1], "native") + unit(2, "pt"), 
              default.units = "native", just = c("left", "bottom"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
  })
  decorate_annotation("Shared clones", {
    grid.text(cs_mtx[co,2], x = seq_along(cs) + 1/6, y = unit(cs_mtx[co,2], "native") + unit(2, "pt"), 
              default.units = "native", just = c("left", "bottom"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
  })
  decorate_annotation("Clones per cell type", {
    grid.text(-ss_mtx[ro,1], x = unit(ss_mtx[ro,1], "native") - unit(2, "pt"), y = rev(seq_along(ss)) + 1/6, 
              default.units = "native", just = c("right", "center"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 0)
  })
  decorate_annotation("Clones per cell type", {
    grid.text(-ss_mtx[ro,2], x = unit(ss_mtx[ro,2], "native") - unit(2, "pt"), y = rev(seq_along(ss)) - 1/6, 
              default.units = "native", just = c("right", "center"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 0)
  })
  dev.off()
}

{
  ## component of CD4T_cko
  
  CD4_TCR_cko_profile <- CD4_TCR_profile %>% filter(`gex:group` == "cko") %>% 
    dplyr::select(CellID, celltype1 = `gex:fine_anno`, clone_id = `airr:clone_id`, clone_id_size = `airr:clone_id_size`)
  CD4_TCR_cko_profile$celltype1 <- factor(CD4_TCR_cko_profile$celltype1, levels = paste0("CD4T_C", 0:5), labels = paste0("C", 1:6))
  
  CD4_TCR_cko_profile <- lapply(unique(CD4_TCR_cko_profile$clone_id), function(i) {
    clones <- CD4_TCR_cko_profile %>% filter(clone_id == i)
    mode <- clones$celltype1 %>% unique() %>% sort() %>% paste(collapse = "-")
    clones$shared_mode <- mode
    clones$code <- mode2code(mode)
    return(clones)
  }) %>% bind_rows()
  
  CD4_TCR_cko_list <- lapply(paste0("C", 1:6), function(i) {
    CD4_TCR_cko_profile %>% filter(celltype1 == i) %>% pull(clone_id)
  })
  names(CD4_TCR_cko_list) <- paste0("C", 1:6)
  
  m_cko = make_comb_mat(CD4_TCR_cko_list, top_n_sets = 6)
  ss <- set_size(m_cko)
  cs <- comb_size(m_cko)
  ss_df <- as.data.frame(table(CD4_TCR_cko_profile$celltype1)) %>% column_to_rownames(colnames(.)[1])
  ss_vec <- setNames(as.numeric(ss_df[names(ss),1]), names(ss))
  cs_df <- as.data.frame(table(CD4_TCR_cko_profile$code)) %>% column_to_rownames(colnames(.)[1])
  cs_vec <- setNames(as.numeric(cs_df[names(cs),1]), names(cs))
  ss_mtx <- matrix(ncol = 2, c(ss_vec,ss)) * -1
  rownames(ss_mtx) <- names(ss)
  cs_mtx <- matrix(ncol = 2, c(cs_vec,cs))
  rownames(cs_mtx) <- names(cs)
  ht <- UpSet(m_cko,
              set_order = 1:6,
              comb_order = order(comb_degree(m_cko), -cs),
              top_annotation = HeatmapAnnotation(
                "Shared clones" = anno_barplot(cs_mtx,
                                               ylim = c(0, max(cs) * 1.1),
                                               border = F,
                                               gp = gpar(fill = c("black", "white"), col = c("black", "black")),
                                               height = unit(4, "cm"),
                                               beside = T, attach = T),
                annotation_name_side = "left",
                annotation_name_rot = 90
              ),
              left_annotation = rowAnnotation(
                "Clones per cell type" = anno_barplot(ss_mtx,
                                                      baseline = 0,
                                                      axis_param = list(
                                                        at = c(0, -500, -1000),
                                                        labels = c(0, 500, 1000),
                                                        labels_rot = 0),
                                                      border = F,
                                                      gp = gpar(fill = c("black", "white"), col = c("black", "black")),
                                                      width = unit(4, "cm"),
                                                      beside = T, attach = T),
                set_name = anno_text(set_name(m_cko),
                                     location = 0.5,
                                     just = "center",
                                     width = max_text_width(set_name(m_cko)) + unit(4, "mm"))
              ),
              right_annotation = NULL,
              show_row_names = F,
              row_names_side = "right",
              width = length(cs)*unit(5, "mm"), height = length(ss)*unit(5, "mm"))
  pdf(file = "results/6_TCR_analysis/CD4_TCR_cko_clone_profile.pdf", width = 15, height = 5)
  ht <- draw(ht)
  co <- column_order(ht)
  ro <- row_order(ht)
  decorate_annotation("Shared clones", {
    grid.text(cs_mtx[co,1], x = seq_along(cs) - 1/6, y = unit(cs_mtx[co,1], "native") + unit(2, "pt"), 
              default.units = "native", just = c("left", "bottom"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
  })
  decorate_annotation("Shared clones", {
    grid.text(cs_mtx[co,2], x = seq_along(cs) + 1/6, y = unit(cs_mtx[co,2], "native") + unit(2, "pt"), 
              default.units = "native", just = c("left", "bottom"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 45)
  })
  decorate_annotation("Clones per cell type", {
    grid.text(-ss_mtx[ro,1], x = unit(ss_mtx[ro,1], "native") - unit(2, "pt"), y = rev(seq_along(ss)) + 1/6, 
              default.units = "native", just = c("right", "center"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 0)
  })
  decorate_annotation("Clones per cell type", {
    grid.text(-ss_mtx[ro,2], x = unit(ss_mtx[ro,2], "native") - unit(2, "pt"), y = rev(seq_along(ss)) - 1/6, 
              default.units = "native", just = c("right", "center"), 
              gp = gpar(fontsize = 6, col = "#404040"), rot = 0)
  })
  dev.off()
}

C3_ctr_profile <- CD4_TCR_ctr_profile %>% 
  filter(celltype1 == "C3") %>%  
  mutate(component = ifelse(shared_mode == "C3", "Single clone",
                            ifelse(str_detect(shared_mode, "C4|C5"), "Shared Treg", "Shared other"))) %>% 
  mutate(component = ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C4)(?=.*C5)"), "Shared C3|C4|C5",
                            ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C4)(?!.*C5)"), "Shared C3|C4", 
                                   ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?!.*C4)(?=.*C5)"), "Shared C3|C5", component)))) %>% 
  mutate(group = "ctr")

C3_cko_profile <- CD4_TCR_cko_profile %>% 
  filter(celltype1 == "C3") %>% 
  mutate(component = ifelse(shared_mode == "C3", "Single clone",
                            ifelse(str_detect(shared_mode, "C4|C5"), "Shared Treg", "Shared other"))) %>% 
  mutate(component = ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C4)(?=.*C5)"), "Shared C3|C4|C5",
                            ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C4)(?!.*C5)"), "Shared C3|C4", 
                                   ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?!.*C4)(?=.*C5)"), "Shared C3|C5", component)))) %>% 
  mutate(group = "cko")

C4_ctr_profile <- CD4_TCR_ctr_profile %>% 
  filter(celltype1 == "C4") %>%  
  mutate(component = ifelse(shared_mode == "C4", "Single clone",
                            ifelse(str_detect(shared_mode, "C3|C5"), "Shared Treg", "Shared other"))) %>% 
  mutate(component = ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C3)(?=.*C5)"), "Shared C3|C4|C5",
                            ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C3)(?!.*C5)"), "Shared C3|C4", 
                                   ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?!.*C3)(?=.*C5)"), "Shared C4|C5", component)))) %>% 
  mutate(group = "ctr")

C4_cko_profile <- CD4_TCR_cko_profile %>% 
  filter(celltype1 == "C4") %>% 
  mutate(component = ifelse(shared_mode == "C4", "Single clone",
                            ifelse(str_detect(shared_mode, "C3|C5"), "Shared Treg", "Shared other"))) %>% 
  mutate(component = ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C3)(?=.*C5)"), "Shared C3|C4|C5",
                            ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C3)(?!.*C5)"), "Shared C3|C4", 
                                   ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?!.*C3)(?=.*C5)"), "Shared C4|C5", component)))) %>% 
  mutate(group = "cko")

C5_ctr_profile <- CD4_TCR_ctr_profile %>% 
  filter(celltype1 == "C5") %>%  
  mutate(component = ifelse(shared_mode == "C5", "Single clone",
                            ifelse(str_detect(shared_mode, "C3|C4"), "Shared Treg", "Shared other"))) %>% 
  mutate(component = ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C3)(?=.*C3)"), "Shared C3|C4|C5",
                            ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C3)(?!.*C4)"), "Shared C3|C5", 
                                   ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?!.*C3)(?=.*C4)"), "Shared C4|C5", component)))) %>% 
  mutate(group = "ctr")

C5_cko_profile <- CD4_TCR_cko_profile %>% 
  filter(celltype1 == "C5") %>% 
  mutate(component = ifelse(shared_mode == "C5", "Single clone",
                            ifelse(str_detect(shared_mode, "C3|C4"), "Shared Treg", "Shared other"))) %>% 
  mutate(component = ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C3)(?=.*C4)"), "Shared C3|C4|C5",
                            ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?=.*C3)(?!.*C4)"), "Shared C3|C5", 
                                   ifelse(component == "Shared Treg" & str_detect(shared_mode, "(?!.*C3)(?=.*C4)"), "Shared C4|C5", component)))) %>% 
  mutate(group = "cko")

C3C4C5_TCR_profile <- do.call(rbind, list(C3_ctr_profile, C3_cko_profile, C4_ctr_profile, C4_cko_profile, C5_ctr_profile, C5_cko_profile))
write_csv(C3C4C5_TCR_profile, file = "results/6_TCR_analysis/CD4_TCR_C3C4C5_profile.csv")

C3C4C5 <- subset(CD4_TCR, subset = celltype1 %in% c("C3", "C4", "C5"))
C3C4C5@meta.data[C3C4C5_TCR_profile$CellID, "celltype2"] <- C3C4C5_TCR_profile$component
C3C4C5$celltype2 <- factor(C3C4C5$celltype2)
C3C4C5$celltype3 <- paste(C3C4C5$celltype1, C3C4C5$celltype2)
C3C4C5$celltype3 <- factor(C3C4C5$celltype3)
saveRDS(C3C4C5, file = "results/6_TCR_analysis/CD4_TCR_C3C4C5.rds")

C3C4C5_ctr <- subset(C3C4C5, subset = group == "ctr" & celltype3 %in% c("C3 Shared C3|C4", "C3 Shared C3|C4|C5", "C3 Shared C3|C5",
                                                                        "C4 Shared C3|C4", "C5 Shared C3|C4|C5", "C5 Shared C3|C5"))
C3C4C5_ctr$celltype3 <- factor(C3C4C5_ctr$celltype3)
levels(C3C4C5_ctr$celltype3) <- c("C3 Shared C3|C4", "C3 Shared C3|C5", "C3 Shared C3|C5", "C4 Shared C3|C4", "C5 Shared C3|C5")

FeatureStatPlot_v5(C3C4C5_ctr, stat.by = c("Foxp3", "Il2ra", "Cst7", "Nkg7", "Gzma", "Prf1", "Tbx21", "Ccl4", "Ccl3", "Ccl5"), 
                   group.by = "celltype3", comparisons = list(c("C3 Shared C3|C4", "C4 Shared C3|C4"),
                                                              c("C3 Shared C3|C5", "C5 Shared C3|C5")))

C3C4C5_cko <- subset(C3C4C5, subset = group == "cko" & celltype3 %in% c("C3 Shared C3|C4", "C3 Shared C3|C4|C5", "C3 Shared C3|C5",
                                                                        "C4 Shared C3|C4", "C5 Shared C3|C4|C5", "C5 Shared C3|C5"))
C3C4C5_cko$celltype3 <- factor(C3C4C5_cko$celltype3)
levels(C3C4C5_cko$celltype3) <- c("C3 Shared C3|C4", "C3 Shared C3|C5", "C3 Shared C3|C5", "C4 Shared C3|C4", "C5 Shared C3|C5", "C5 Shared C3|C5")

FeatureStatPlot_v5(C3C4C5_cko, stat.by = c("Foxp3", "Il2ra", "Cst7", "Nkg7", "Gzma", "Prf1", "Tbx21", "Ccl4", "Ccl3", "Ccl5"), 
                   group.by = "celltype3", comparisons = list(c("C3 Shared C3|C4", "C4 Shared C3|C4"),
                                                              c("C3 Shared C3|C5", "C5 Shared C3|C5")))

# 计算exTreg score
library(GSEABase)

exTreg_sig <- list(exTreg = c("Cst7", "Nkg7", "Gzma", "Prf1", "Tbx21", "Ccl4"))
C3C4C5_sub <- subset(C3C4C5, subset = celltype3 %in% c("C3 Shared C3|C4", "C3 Shared C3|C4|C5", "C3 Shared C3|C5",
                                                       "C4 Shared C3|C4", "C5 Shared C3|C4|C5", "C5 Shared C3|C5"))
C3C4C5_sub$celltype3 <- factor(C3C4C5_sub$celltype3)
levels(C3C4C5_sub$celltype3) <- c("C3 Shared C3|C4", "C3 Shared C3|C5", "C3 Shared C3|C5",
                                  "C4 Shared C3|C4", "C5 Shared C3|C5", "C5 Shared C3|C5")
C3C4C5 <- AddModuleScore_UCell(C3C4C5, features = exTreg_sig, name = NULL, ncores = 12)
score_c <-  C3C4C5@meta.data %>% dplyr::select(c("group", "celltype3", names(exTreg_sig))) %>% rownames_to_column("CellID")
write_csv(score_c, file = "results/6_TCR_analysis/CD4_TCR_C3C4C5_score_c.csv")

C3C4C5_sub1 <- subset(C3C4C5_sub, subset = celltype3 %in% c("C3 Shared C3|C4", "C3 Shared C3|C5"))
C3C4C5_sub1$celltype3 <- factor(C3C4C5_sub1$celltype3)
FeatureStatPlot_v5(C3C4C5_sub1, stat.by = c("Foxp3", "Il2ra"), group.by = "celltype3", split.by = "group", comparisons = T)
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_C3shared_ctrvscko.pdf", width = 6, height = 3)
C3C4C5_sub2 <- subset(C3C4C5_sub, subset = celltype3 %in% c("C3 Shared C3|C4", "C4 Shared C3|C4"))
C3C4C5_sub2$celltype3 <- factor(C3C4C5_sub2$celltype3)
FeatureStatPlot_v5(C3C4C5_sub2, stat.by = c("Foxp3", "Il2ra"), group.by = "group", split.by = "celltype3", comparisons = T)
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_shared_C3vsC4.pdf", width = 8, height = 3)
C3C4C5_sub3 <- subset(C3C4C5_sub, subset = celltype3 %in% c("C3 Shared C3|C5", "C5 Shared C3|C5"))
C3C4C5_sub3$celltype3 <- factor(C3C4C5_sub3$celltype3)
FeatureStatPlot_v5(C3C4C5_sub3, stat.by = c("Foxp3", "Il2ra"), group.by = "group", split.by = "celltype3", comparisons = T)
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_shared_C3vsC5.pdf", width = 8, height = 3)

C3 <- subset(C3C4C5, subset = celltype1 == "C3")
C3$celltype3 <- factor(C3$celltype3)
FeatureStatPlot_v5(C3C4C5, stat.by = c("Ifng"), group.by = "group", split.by = "celltype3", comparisons = T)

# 计算C3中exTreg及Th1-like的统计量
df <- GetAssayData(C3C4C5, layer = "data") %>% .[c("Cd40lg", "Tbx21", "Ifng", "Foxp3", "Il2ra"),] %>% 
  as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("CellID") %>% 
  left_join(score_c, ., by = "CellID")
df <- df %>% filter(str_detect(celltype3, "^C3|^C4 Shared C3"))
df$celltype3 <- factor(df$celltype3)
levels(df$celltype3) <- c("C3 Shared C3|C4", "C3 Shared C3|C5", "C3 Shared C3|C5", "C3 Others", "C3 Others",
                          "C4 Shared C3|C4", "C4 Shared C3|C4")
df$celltype3 <- factor(df$celltype3, levels = c("C3 Others", "C3 Shared C3|C5", "C3 Shared C3|C4", "C4 Shared C3|C4"))

# Kruskal-Wallis test
kruskal_results <- df %>%
  group_by(group) %>%
  summarise(p.value = kruskal.test(Ifng ~ celltype3)$p.value)
# Dunn's 检验
dunn_results <- df %>%
  group_by(group) %>%
  do(dunnTest(Ifng ~ celltype3, data = ., method = "bonferroni")$res)

library(ggstatsplot)
library(tinyfuncr)
library(ggbeeswarm)
ggbetweenstats(data = filter(df, group == "cko"),
               x = celltype3, y = Tbx21,
               plot.type = "violin",
               type = "nonparametric",
               pairwise.display = "none",
               p.adjust.method = "bonferroni",
               results.subtitle = F,
               centrality.plotting = T,
               point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.8), alpha =
                                   0.6, size = 1, stroke = 0, na.rm = TRUE),
               boxplot.args = list(width = 0.01, alpha = 0.2, na.rm = TRUE),
               violin.args = list(width = 1, alpha = 0.2, na.rm = TRUE),
               ggsignif.args = list(textsize = 3, tip_length = 0.01, na.rm = TRUE),
               ggtheme = theme_niwot(),
               package = "ggsci",
               palette = "nrc_npg") +
  theme(axis.text.x = element_blank())+
  scale_color_manual(values = c("#474747", "#3C5488", "#DC0000", "#096D30")) +
  labs(x = "", y = "Tbx21")

ggplot(df, aes(x = celltype3, y = exTreg, color = celltype3)) +
  geom_violin(trim = T, width = 1, color = "black", lwd = 0.3, fill = "transparent", position = position_dodge(0.8), alpha = 0.2, na.rm = T) +
  # geom_boxplot(width = 0.1, color = "black", lwd = 0.3, position = position_dodge(0.8), alpha = 0.2, outlier.shape = NA, na.rm = T) +
  geom_jitter(alpha = 0.4, size = 0.5, width = 0.2) +
  geom_signif(
    comparisons = list(c("C3 Shared C3|C5", "C3 Others"),
                       c("C3 Shared C3|C4", "C3 Shared C3|C5"),
                       c("C3 Shared C3|C4", "C3 Others"),
                       c("C3 Shared C3|C4", "C4 Shared C3|C4")),
    map_signif_level = T,
    step_increase = 0.1,
    tip_length = 0.02,
    size = 0.3,
    textsize = 2
  ) +
  facet_grid(. ~ group) +
  theme_niwot() + 
  theme(legend.position = "right", 
        axis.text.x = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = c("#474747", "#2878B5", "#C82423", "#096D30")) +
  labs(x = "", y = "exTreg")
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_C3shared_exTreg.pdf", width = 6, height = 3)
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_C3shared_Cd40lg.pdf", width = 6, height = 3)
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_C3shared_Tbx21.pdf", width = 6, height = 3)
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_C3shared_Ifng.pdf", width = 6, height = 3)

# 计算C4中Foxp3及Il2ra的统计量
df2 <- GetAssayData(C3C4C5, layer = "data") %>% .[c("Foxp3", "Il2ra"),] %>% 
  as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("CellID") %>% 
  left_join(score_c, ., by = "CellID")
df2 <- df2 %>% filter(str_detect(celltype3, "^C4") | celltype3 == "C3 Shared C3|C4")
df2$celltype3 <- factor(df2$celltype3)
levels(df2$celltype3) <- c("C3 Shared C3|C4", "C4 Shared C3|C4", "C4 Shared C3|C4", "C4 Shared C4|C5", "C4 Others",
                          "C4 Others")
df2$celltype3 <- factor(df2$celltype3, levels = c("C4 Others", "C4 Shared C4|C5", "C4 Shared C3|C4", "C3 Shared C3|C4"))


# Kruskal-Wallis test
kruskal_results <- df2 %>%
  group_by(group) %>%
  summarise(p.value = kruskal.test(Il2ra ~ celltype3)$p.value)
# Dunn's 检验
dunn_results <- df2 %>%
  group_by(group) %>%
  do(dunnTest(Il2ra ~ celltype3, data = ., method = "bonferroni")$res)

ggplot(df2, aes(x = celltype3, y = exTreg, color = celltype3)) +
  geom_violin(trim = T, width = 1, color = "black", lwd = 0.3, fill = "transparent", position = position_dodge(0.8), alpha = 0.2, na.rm = T) +
  # geom_boxplot(width = 0.1, color = "black", lwd = 0.3, position = position_dodge(0.8), alpha = 0.2, outlier.shape = NA, na.rm = T) +
  geom_jitter(alpha = 0.4, size = 0.5, width = 0.2) +
  geom_signif(
    comparisons = list(c("C3 Shared C3|C4", "C4 Shared C3|C4"),
                       c("C3 Shared C3|C4", "C4 Shared C4|C5"),
                       c("C3 Shared C3|C4", "C4 Others"),
                       c("C4 Shared C3|C4", "C4 Others")),
    map_signif_level = T,
    step_increase = 0.1,
    tip_length = 0.02,
    size = 0.3,
    textsize = 2
  ) +
  facet_grid(. ~ group) +
  theme_niwot() + 
  theme(legend.position = "right", 
        axis.text.x = element_blank(),
        strip.background = element_blank()) +
  scale_color_manual(values = c("#474747", "#2878B5", "#096D30", "#C82423")) +
  labs(x = "", y = "exTreg")

ggsave(filename = "results/6_TCR_analysis/CD4_TCR_C4shared_exTreg.pdf", width = 6, height = 3)
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_C4shared_Foxp3.pdf", width = 6, height = 3)
ggsave(filename = "results/6_TCR_analysis/CD4_TCR_C4shared_Il2ra.pdf", width = 6, height = 3)

# C3 component statistics
# 数据准备
group_A <- c(56, 18, 167)  # A组中红球、蓝球、绿球数量
group_B <- c(264, 134, 448)  # B组中红球、蓝球、绿球数量

# 计算总数
total_A <- sum(group_A)
total_B <- sum(group_B)

# 使用卡方检验比较整体差异
chisq.test(rbind(group_A, group_B))

# Shared with Tregs
chisq.test(matrix(c(group_A[1], group_B[1], sum(group_A) - group_A[1], sum(group_B) - group_B[1]), ncol=2))
fisher.test(matrix(c(group_A[1], group_B[1], sum(group_A) - group_A[1], sum(group_B) - group_B[1]), ncol=2))

### 11. CellChat analysis ----

## 数据准备
C3C4C5$celltype4 <- factor(C3C4C5$celltype3)
levels(C3C4C5$celltype4) <- c("C3 Shared C3|C4", "C3 Shared C3|C4", "C3 Shared C3|C5", "C3 Others", "C3 Others",
                              "C4 Shared C3|C4", "C4 Shared C3|C4", "C4 Shared C4|C5", "C4 Others", "C4 Others",
                              "C5 Shared C3|C5", "C5 Shared C3|C5", "C5 Shared C4|C5", "C5 Others", "C5 Others")
C2C3C4 <- subset(CD8T, subset = celltype1 %in% c("C2", "C3", "C4"))
C2C3C4$celltype2 <- factor(C2C3C4$celltype1)
levels(C2C3C4$celltype2) <- paste0("CD8_", c("C2", "C3", "C4"))
DC <- subset(alldata, subset = fine_anno %in% c("DC_C0", "DC_C1", "DC_C2"))
DC$fine_anno <- factor(DC$fine_anno)
cellchat_input <- alldata[,c(colnames(C3C4C5), colnames(C2C3C4), colnames(DC))]
cellchat_input$celltype1 <- ""
cellchat_input@meta.data[colnames(C3C4C5),"celltype1"] <- as.character(C3C4C5@meta.data[colnames(C3C4C5), "celltype4"])
cellchat_input@meta.data[colnames(C2C3C4),"celltype1"] <- as.character(C2C3C4@meta.data[colnames(C2C3C4), "celltype2"])
cellchat_input@meta.data[colnames(DC),"celltype1"] <- as.character(DC@meta.data[colnames(DC), "fine_anno"])
table(cellchat_input$celltype1)
cellchat_input$celltype1 <- factor(cellchat_input$celltype1)
saveRDS(cellchat_input, file = "results/7_CellChat_analysis/cellchat_input.rds")

#### MoMa ----
rm(list = ls()[-which(ls() %in% c("alldata"))])
gc()
MoMa <- subset(alldata, subset = main_anno == "Macro")
MoMa$main_anno <- factor(MoMa$main_anno)
MoMa$fine_anno <- factor(MoMa$fine_anno)
CellDimPlot(MoMa, group.by = "fine_anno", reduction = "umap.scanorama", label = T)
CellStatPlot(MoMa, stat.by = "fine_anno", group.by = "group", plot_type = "trend", label = T)

## anno rename
levels(MoMa$fine_anno) <- paste0("MoMa_C", c("0", "3", "1", "2", "4"))
MoMa$fine_anno <- factor(MoMa$fine_anno, levels = paste0("MoMa_C", c("0", "1", "2", "3", "4")))

## plot UMAP and Proportion
CellDimPlot(MoMa, group.by = "fine_anno", reduction = "umap.scanorama", label = T, palcolor = colors[c(6,5,4,1,2)])
ggsave(filename = "results/4_funtional_analysis/12_MoMa_cluster.pdf", height = 5, width = 5)
CellDimPlot(MoMa, group.by = "fine_anno", reduction = "umap.scanorama", split.by = "group", label = T, palcolor = colors[c(6,5,4,1,2)])
ggsave(filename = "results/4_funtional_analysis/12_MoMa_cluster_by_group.pdf", height = 5, width = 10)
CellStatPlot(MoMa, stat.by = "fine_anno", group.by = "group", plot_type = "trend", label = T, palcolor = colors[c(6,5,4,1,2)])
ggsave(filename = "results/4_funtional_analysis/12_MoMa_anno_proportion.pdf", width = 3.5, height = 5)

## markers
FeatureDimPlot_v5(MoMa, features = c("Itgam", "Csf1r", "Cd68", "Cd86", "Cd80", "Adgre1", "Fcgr4", "Mrc1", "Cd163"), reduction = "umap.scanorama") # Macrophage
FeatureDimPlot_v5(MoMa, features = c("Lair1", "Havcr2", "Lgals9", "Vsir"), reduction = "umap.scanorama") # Macrophage
FeatureDimPlot_v5(MoMa, features = c("Itgam", "Ly6c1", "Ly6g", "Ly6c2", "Ccr2", "Cd14", "Fcgr3", "Adgre1", "Arg1", "Ms4a7", "Cst3"), 
            reduction = "umap.scanorama") # 区分macrophage及monocyte
FeatureDimPlot_v5(MoMa, features = c("Nos2", "Cd1d1", "Hopx", "Nr4a2", "Bhlhe40"), reduction = "umap.scanorama", split.by = "group", ncol = 2, color_limits = c(0,3))
ggsave(filename = "results/4_funtional_analysis/12_MoMa_gene1.pdf", width = 6, height = 10)
FeatureDimPlot_v5(MoMa, features = c("Il1rn", "Il1a", "Il1b", "Il12a", "Il12b", "Il12rb1", "Il12rb2", "Il23a", "Il23r", "Il27", "Ebi3", "Il27ra", "Il6st"), 
                  reduction = "umap.scanorama", split.by = "group")
ggsave(filename = "results/4_funtional_analysis/12_MoMa_gene2.pdf", width = 18, height = 10)
 
# https://www.nature.com/articles/s41467-021-21407-w
FeatureDimPlot_v5(MoMa, features = c("Ly6c2", "Plac8", "Ccr2", "Tgfbi", "Chil3", "Ace", "Ifitm2", "Ifitm3", "S100a6"), reduction = "umap.scanorama") # monocyte signature
FeatureDimPlot_v5(MoMa, features = c("Cxcl9", "Spp1"), reduction = "umap.scanorama")

FeatureDimPlot_v5(MoMa, features = c("Nos2", "Mrc1", "Ccr2", "Cx3cr1", "Irf5", "Cd274", "Cd1d1"), reduction = "umap.scanorama", split.by = "group")
FeatureDimPlot_v5(MoMa, features = c("Vegfa", "Ptgs2", "H2-Eb1"), reduction = "umap.scanorama", split.by = "group")
FeatureDimPlot_v5(MoMa, features = c("Cd80", "Cd86", "Cd40"), reduction = "umap.scanorama", split.by = "group")

## GSEA analysis
C0C4 <- subset(MoMa, subset = fine_anno %in% c("MoMa_C0", "MoMa_C4"))
C0C4$fine_anno <- factor(C0C4$fine_anno)
Idents(C0C4) <- "fine_anno"
C0vsC4 <- FindMarkers(C0C4, ident.1 = "MoMa_C0", ident.2 = "MoMa_C4", logfc.threshold = 0)
C0vsC4 <- C0vsC4 %>% rownames_to_column("ID")
write_csv(C0vsC4, file = "results/4_funtional_analysis/12_MoMa_C0vsC4.csv")

C0 <- subset(MoMa, subset = fine_anno %in% c("MoMa_C0"))
C0$fine_anno <- factor(C0$fine_anno)
Idents(C0) <- "group"
C0_cKOvsctr <- FindMarkers(C0, ident.1 = "cko", ident.2 = "ctr", logfc.threshold = 0)
C0_cKOvsctr <- C0_cKOvsctr %>% rownames_to_column("ID")
write_csv(C0_cKOvsctr, file = "results/4_funtional_analysis/12_MoMa_C0_cKOvsctr.csv")

C4 <- subset(MoMa, subset = fine_anno %in% c("MoMa_C4"))
C4$fine_anno <- factor(C4$fine_anno)
Idents(C4) <- "group"
C4_cKOvsctr <- FindMarkers(C4, ident.1 = "cko", ident.2 = "ctr", logfc.threshold = 0)
C4_cKOvsctr <- C4_cKOvsctr %>% rownames_to_column("ID")
write_csv(C4_cKOvsctr, file = "results/4_funtional_analysis/12_MoMa_C4_cKOvsctr.csv")

## M1,M2 signature
library(GSEABase)
## homologous gene conversion
Homologous_gene <- read_tsv(file = "~/Bioinformatics/Programs/2023/02_yanglichao_scTransSpecies/mart_export.txt") %>% dplyr::select(1:9)

immu_gmt <- getGmt("~/Bioinformatics/Programs/2022/8_public_database/data/DIY_gmt/hsa_Immune_response.gmt")
genesets <- lapply(names(immu_gmt), function(x) {
  gene <- immu_gmt[[x]]@geneIds
  tmp <- Homologous_gene %>% filter(`Gene name` %in% gene) %>% pull(`Mouse gene name`) %>% 
    na.omit() %>% unique()
  tmp <- intersect(tmp, rownames(MoMa))
  return(tmp)
})
names(genesets) <- names(immu_gmt)
genesets <- genesets[c(17:18,30:33,49:52,69:76,78:86,101:111)]

library(UCell)
MoMa <- AddModuleScore_UCell(MoMa, features = genesets, name = NULL, ncores = 12)
score_c <- MoMa@meta.data %>% dplyr::select(c("group", "fine_anno", names(genesets)))
write_csv(score_c, file = "results/4_funtional_analysis/12_Macro_signature_Ucell.csv")

# statistics analysis
library(ggsignif)
signature_meta <- data.frame(Signatures = names(genesets), Labels = paste0("S", 1:length(genesets)))
data <- score_c %>% dplyr::select(-1) %>% filter(fine_anno %in% c("MoMa_C0", "MoMa_C4"))1
colnames(data) <- c("fine_anno", signature_meta$Labels)
data <- data %>% filter(S29 != 0) %>% mutate(M1_M2_ratio = S28 / S29)

# data transforming
data <- data %>% 
  pivot_longer(cols = -fine_anno, names_to = "Signature", values_to = "Scores")
data$fine_anno <- factor(data$fine_anno, levels = c("MoMa_C0", "MoMa_C4"))
data$Signature <- factor(data$Signature, levels = paste0("S", 1:length(genesets)))

# Wilcoxon test
wilcoxon_results <- data %>%
  group_by(Signature) %>%
  summarise(p.value = wilcox.test(Scores ~ fine_anno)$p.value) %>% left_join(signature_meta, by = c("Signature" = "Labels"))

p <- lapply(c(paste0("S",1:38), "M1_M2_ratio"), function(i) {
  df1 <- data %>% filter(Signature == i)
  ggplot(df1, aes(x = Signature, y = Scores, fill = fine_anno)) +
    geom_violin(trim = F, width = 0.8, position = position_dodge(0.8), alpha = 0.5) +
    geom_boxplot(width = 0.08, position = position_dodge(0.8), outlier.shape = NA) +
    theme_niwot() + theme(legend.position = "top") +
    scale_fill_manual(values = c("MoMa_C0" = "#DC0000", "MoMa_C4" = "#3C5488")) +
    labs(x = "", y = "Scores")
})
pdf(file = "results/4_funtional_analysis/12_signature_all.pdf", width = 4, height = 4)
print(p)
dev.off()
