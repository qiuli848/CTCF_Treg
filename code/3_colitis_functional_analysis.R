rm(list = ls())
gc()

library(reticulate)
use_condaenv("scrna", required = TRUE)

library(tidyverse)
library(Seurat)
library(cowplot)
library(SCP)
library(ggsci)
library(colorspace)
library(RColorBrewer)
source("/home/data/qiuli/Bioinformatics/code2source/SCP_adjustment_for_v5.R", echo=TRUE)

alldata <- readRDS(file = "results/3_anno.rds")
CellDimPlot(alldata, group.by = "main_anno", reduction = "umap.scanorama", label = T)

#### CD4T ----
rm(list = ls()[-which(ls() %in% c("alldata"))])
gc()
CD4T <- subset(alldata, subset = main_anno == "CD4T" | fine_anno == "CD4T_Cycling")
CD4T$main_anno <- factor(CD4T$main_anno)
CD4T$fine_anno <- factor(CD4T$fine_anno)
# 重新降维
CD4T <- CD4T %>% 
  # RunPCA(npcs = 30, verbose = F) %>% 
  RunUMAP(dims = 1:30, reduction = "pca")
CellDimPlot(CD4T, group.by = "fine_anno", label = T)
CD4T$main_anno <- factor(CD4T$main_anno)
CD4T$fine_anno <- factor(CD4T$fine_anno)
CD4T$group <- factor(CD4T$group, levels = c("CPI_normal", "CPI_nocolitis", "CPIcolitis"))
CD4T$celltype1 <- factor(CD4T$fine_anno)
levels(CD4T$celltype1) <- paste0("C", 1:9)
saveRDS(CD4T, file = "results/CD4T.rds")

CD4T@meta.data <- CD4T@meta.data %>%
  mutate(fine_anno = case_when(
    fine_anno == "CD4T_C1_Treg_FOXP3" ~ "CD4T_C2_Treg_FOXP3",
    fine_anno == "CD4T_C2_Treg_CCL5" ~ "CD4T_C3_Treg_like_CCL5",
    fine_anno == "CD4T_C3_Treg_TOX2" ~ "CD4T_C1_Tex_TOX2",
    TRUE ~ fine_anno
  ))
CD4T$fine_anno <- factor(CD4T$fine_anno)
CD4T$celltype1 <- factor(CD4T$fine_anno)
levels(CD4T$celltype1) <- paste0("C", 1:9)
saveRDS(CD4T, file = "results/CD4T.rds")

Idents(CD4T) <- "celltype1"
CD4_markers <- FindAllMarkers(CD4T)
write_csv(CD4_markers, file = "results/6_functional_analysis/6_CD4T_markers.csv")

c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")

### 6. UMAP
colors <- brewer.pal(12, "Paired")
CellDimPlot(CD4T, group.by = "celltype1", label = T, palcolor = colors[c(1,2,6,4,5,3,9,7,8)])
ggsave(filename = "results/6_functional_analysis/6_CD4T_cluster.pdf", height = 5, width = 5)
CellDimPlot(CD4T, group.by = "celltype1", split.by = "group", palcolor = colors[c(1,2,6,4,5,3,9,7,8)])
ggsave(filename = "results/6_functional_analysis/6_CD4T_cluster_by_group.pdf", height = 5, width = 15)

### 6. cell proportion
CellStatPlot(CD4T, stat.by = "fine_anno", group.by = "group", plot_type = "trend", label = T, palcolor = colors[c(1,2,6,4,5,3,9,7,8)])
ggsave(filename = "results/6_functional_analysis/6_CD4T_proportion.pdf", width = 6, height = 6)

### 6. marker gene umap
FeatureDimPlot_v5(CD4T, features = c("FOXP3", "IL2RA", "CST7", "NKG7", "GZMA", "PRF1", "TBX21", "CCL4", "CCL3", "CCL5", "CD40LG", "TNFSF9"))
ggsave(filename = "results/6_functional_analysis/6_exTreg_gene_umap.pdf", height = 9, width = 12)

FeatureDimPlot_v5(CD4T, features = c("CD40LG", "TNFSF9"), split.by = "group")
ggsave(filename = "results/6_functional_analysis/6_exTreg_CD40L_umap_by_group.pdf", height = 7, width = 10)

FeatureDimPlot_v5(alldata, features = c("CXCR5", "CXCL13"), split.by = "group")
FeatureDimPlot_v5(CD4T, features = c("IKZF2", "GATA3", "RORC", "MAF"), split.by = "group", ncol = 3)

### 6. Differential expression analysis
Idents(CD4T) <- "celltype1"
CD4T_C3vsC7 <- FindMarkers(CD4T, ident.1 = "C3", ident.2 = "C7", logfc.threshold = 0) %>% rownames_to_column("gene")
write_csv(CD4T_C3vsC7, file = "results/6_functional_analysis/6_CD4T_exTregvsTh1.csv")
CD4T_C3vsC2 <- FindMarkers(CD4T, ident.1 = "C3", ident.2 = "C2", logfc.threshold = 0) %>% rownames_to_column("gene")
write_csv(CD4T_C3vsC2, file = "results/6_functional_analysis/6_CD4T_exTregvsTreg.csv")

#### CD8T
rm(list = ls()[-which(ls() %in% c("alldata"))])
gc()
CD8T <- subset(alldata, subset = main_anno == "CD8T" | main_anno == "MAIT" | fine_anno == "CD8T_Cycling")
CD8T$main_anno <- factor(CD8T$main_anno)
CD8T$fine_anno <- factor(CD8T$fine_anno)
CellDimPlot(CD8T, group.by = "fine_anno", reduction = "umap.scanorama", label = T)

CD8T <- CD8T %>% 
  # RunPCA(npcs = 30, verbose = F) %>%
  RunUMAP(dims = 1:30, reduction = "pca")
CellDimPlot(CD8T, group.by = "fine_anno", label = T)
CD8T$main_anno <- factor(CD8T$main_anno)
CD8T$fine_anno <- factor(CD8T$fine_anno)
CD8T$group <- factor(CD8T$group, levels = c("CPI_normal", "CPI_nocolitis", "CPIcolitis"))
CD8T$celltype1 <- factor(CD8T$fine_anno)
levels(CD8T$celltype1) <- paste0("C", 0:6)
Idents(CD8T) <- "celltype1"

### 6. UMAP
colors <- brewer.pal(12, "Paired")
CellDimPlot(CD8T, group.by = "celltype1", label = T, palcolor = colors)
ggsave(filename = "results/6_functional_analysis/1_CD8T_cluster.pdf", height = 5, width = 5)
CellDimPlot(CD8T, group.by = "celltype1", split.by = "group", palcolor = colors)
ggsave(filename = "results/6_functional_analysis/1_CD8T_cluster_by_group.pdf", height = 5, width = 15)

### 6. cell proportion
CellStatPlot(CD8T, stat.by = "fine_anno", group.by = "group", plot_type = "trend", label = T, palcolor = colors)
ggsave(filename = "results/6_functional_analysis/1_CD8T_proportion.pdf", width = 6, height = 6)

#### myeloid ----
rm(list = ls()[-which(ls() %in% c("alldata"))])
gc()

myeloid <- subset(alldata, subset = main_anno == "Myeloid" & fine_anno != "pDC")2020310771

myeloid$main_anno <- factor(myeloid$main_anno)
myeloid$fine_anno <- factor(myeloid$fine_anno)
myeloid$group <- factor(myeloid$group, levels = c("CPI_normal", "CPI_nocolitis", "CPIcolitis"))
CellDimPlot(myeloid, group.by = "fine_anno", reduction = "umap.scanorama", label = T)

myeloid <- myeloid %>% 
  # RunPCA(npcs = 30, verbose = F) %>%
  RunUMAP(dims = 1:30, reduction = "pca") %>% 
  FindNeighbors(reduction = "pca", dims = 1:30) %>% 
  FindClusters(graph.name = "RNA_snn", resolution = 0.3, algorithm = 1)
CellDimPlot(myeloid, group.by = "RNA_snn_res.0.3", label = T)
CellDimPlot(myeloid, group.by = "fine_anno", label = T)

myeloid$celltype1 <- factor(myeloid$RNA_snn_res.0.3)
levels(myeloid$celltype1) <- c("Mono", "Macro1", "Macro2", "DC", "Macro3")
myeloid$celltype1 <- factor(myeloid$celltype1, levels = c("DC", "Mono", "Macro1", "Macro2", "Macro3"))
Idents(myeloid) <- "celltype1"
CellDimPlot(myeloid, group.by = "celltype1", label = T, palcolor = colors)
ggsave(filename = "results/6_functional_analysis/12_myeloid_cluster.pdf", height = 5, width = 5)

FeatureStatPlot_v5(myeloid, stat.by = c("CD40", "TNFRSF9"), group.by = "celltype1", split.by = "group", comparisons = T)
ggsave(filename = "results/6_functional_analysis/12_myeloid_CD40_group.pdf", height = 3, width = 12)

library(SeuratObject)
library(FSA) # Dunn's test
df <- as.matrix(GetAssayData(myeloid, layer = "data")["CD40",]) %>% 
  as.data.frame() %>% rownames_to_column("CellID") %>% rename(CD40 = V1) %>% 
  left_join(myeloid@meta.data %>% filter(celltype1 %in% c("DC", "Mono")) %>% select(celltype1, group) %>% rownames_to_column("CellID"),
            ., by = "CellID")

# Kruskal-Wallis test
kruskal_results <- df %>%
  group_by(celltype1) %>%
  summarise(p.value = kruskal.test(CD40 ~ group)$p.value)
# Dunn's 检验
dunn_results <- df %>%
  group_by(celltype1) %>%
  do(dunnTest(CD40 ~ group, data = ., method = "bonferroni")$res)

### 12. cell proportion
CellStatPlot(myeloid, stat.by = "celltype1", group.by = "group", plot_type = "trend", label = T, palcolor = colors)
ggsave(filename = "results/6_functional_analysis/12_myeloid_proportion.pdf", width = 5, height = 5)

markers <- FindAllMarkers(myeloid)

FeatureDimPlot_v5(myeloid, features = c("CD80", "CD86", "CD40", "CD83"), split.by = "group", ncol = 3)
FeatureDimPlot_v5(myeloid, features = c("CCR2", "CX3CR1"), split.by = "group", ncol = 3)
levels(myeloid$celltype1) <- c("Mono_FCN1", "Macro1_ACP5", "Macro2_RNASE1", "DC_FLT3", "Macro3_phag.")

mono <- subset(myeloid, subset = celltype1 == "Mono_FCN1")
Idents(mono) <- "group"
mono_colitisvsno <- FindMarkers(mono, ident.1 = "CPIcolitis", ident.2 = "CPI_nocolitis", logfc.threshold = 0) %>% 
  rownames_to_column("gene")
write_csv(mono_colitisvsno, file = "results/6_functional_analysis/12_mono_colitisvsno.csv")

DC <- subset(myeloid, subset = celltype1 == "DC_FLT3")
Idents(DC) <- "group"
DC_colitisvsno <- FindMarkers(DC, ident.1 = "CPIcolitis", ident.2 = c("CPI_nocolitis", "CPI_normal"), logfc.threshold = 0) %>% 
  rownames_to_column("gene")
write_csv(DC_colitisvsno, file = "results/6_functional_analysis/12_DC_colitisvsother.csv")

#### TCR analysis
## convert seuratobj to anndata
adata <- srt_to_adata_v5(srt = alldata, 
                         meta_to_use = c("orig.ident", "nCount_RNA", "nFeature_RNA", "sample", "group", "main_anno", "fine_anno"),
                         reduction_to_use = c("scanorama", "umap.scanorama"))
## Save data
adata$write_h5ad("results/5_scTCR/alldata_to_anndata_raw.h5ad")

bar_data <- read_csv(file = "results/5_scTCR/TCR_profile.csv")
bar_data <- bar_data %>% filter(`gex:main_anno` %in% c("CD4T", "Cycling") & `gex:fine_anno` != "CD8T_Cycling") %>% 
  dplyr::select(group = `gex:group`, celltype = `gex:fine_anno`, clonal_expansion = `airr:clonal_expansion`)
write_csv(bar_data, file = "results/5_scTCR/bar_data.csv")

