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
library(ggpubr)
setwd("/home/data/qiuli/Bioinformatics/Programs/2025/02_public_PD1_scRNA_scTCR/analysis")

paired_color <- c(
  "#A6CEE3","#1F78B4","#FB9A99","#E31A1C","#B2DF8A","#33A02C",
  "#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928",
  "#B3E2CD","#1B9E77","#F4CAE4","#E7298A"
)

CD4 <- readRDS(file = "3_anno_CD4.rds")
CellDimPlot(CD4, group.by = "celltype1", reduction = "umap.harmony", palcolor = paired_color, raster = T)
ggsave(filename = "results/3_analysis/01_CD4_UMAP.pdf", height = 5)
FeatureDimPlot(CD4, features = c("FOXP3", "IL2RA", "IKZF2", "NOTCH3"), reduction = "umap.harmony", raster = T)

library(UCell)
library(FSA)
library(GSEABase)

Hallmarker <- getGmt("/home/data/qiuli/Bioinformatics/data_public/msigdb_v2023.1.Hs_GMTs/h.all.v2023.1.Hs.symbols.gmt")
GO_BP <- getGmt("/home/data/qiuli/Bioinformatics/data_public/msigdb_v2023.1.Hs_GMTs/c5.go.bp.v2023.1.Hs.symbols.gmt")

genesets <- list(
  IFN_sig = Hallmarker[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]@geneIds,
  TregInst = c("FOXP3", "CTLA4", "TNFRSF9", "TIGIT", "IKZF2", "FGL2"),
  Th1_sig = GO_BP[["GOBP_T_HELPER_1_CELL_DIFFERENTIATION"]]@geneIds,
  Th2_sig = GO_BP[["GOBP_T_HELPER_2_CELL_DIFFERENTIATION"]]@geneIds,
  Type1_ck = GO_BP[["GOBP_REGULATION_OF_T_HELPER_1_TYPE_IMMUNE_RESPONSE"]]@geneIds,
  Treg_sig = GO_BP[["GOBP_REGULATORY_T_CELL_DIFFERENTIATION"]]@geneIds
)
CD4 <- AddModuleScore_UCell(CD4, features = genesets, name = NULL, ncores = 12)
score_c <-  CD4@meta.data %>% dplyr::select(c("sample_id", "celltype2", names(genesets))) %>% rownames_to_column("CellID")
write_csv(score_c, file = "results/3_analysis/CD4_score_c_genesets.csv")

Treg <- subset(CD4, subset = celltype2 %in% c("C1", "C2", "C3", "C4", "C5"))
Tcon <- subset(CD4, subset = celltype2 %in% c("C6", "C7", "C8", "C9", "C10", "C11" ,"C12", "C13", "C14", "C15"))
Treg$celltype2 <- factor(Treg$celltype2, levels = c("C1", "C2", "C3", "C4", "C5"))


### 1.细胞成分比较 ----

# 提取meta.data信息
meta <- CD4@meta.data %>%
  filter(dataset == "GSE176021") %>%
  filter(celltype2 %in% c("C1", "C2", "C3", "C4", "C5", "C8", "C9", "C11", "C12", "C13")) %>%
  dplyr::select(sample_id, response, treatment, celltype2)

# 计算 sample-level 细胞比例
fraction <- meta %>%
  filter(response %in% c("NR", "R")) %>%
  group_by(sample_id, response, treatment, celltype2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample_id, response, treatment) %>%
  mutate(Fraction = n / sum(n)) %>%
  ungroup()
fraction$response <- factor(fraction$response, levels = c("NR", "R"))
fraction$treatment <- factor(fraction$treatment, levels = c("pre", "post"))

# 作图风格设置
theme_nature <- theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.line = element_line(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

# 作图
p1 <- ggplot(
  fraction,
  aes(
    x = celltype2,
    y = Fraction,
    fill = response
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65,
    alpha = 0.6,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color = response),
    size = 1,
    alpha = 0.9,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)
  ) +
  stat_compare_means(
    aes(group = response),
    method = "wilcox.test",
    label = "p.format",
    size = 3
  ) +
  facet_wrap(~ treatment, ncol = 1) +
  scale_fill_manual(values = c("NR" = "#B2182B", "R" = "#2166AC")) +
  scale_color_manual(values = c("NR" = "#B2182B", "R" = "#2166AC")) +
  labs(
    x = NULL,
    y = "Fraction of CD4 T cells"
  ) +
  theme_nature
print(p1)

p2 <- ggplot(
  fraction,
  aes(
    x = celltype2,
    y = Fraction,
    fill = treatment
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.65,
    alpha = 0.6,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color = treatment),
    size = 1,
    alpha = 0.9,
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    )
  ) +
  stat_compare_means(
    aes(group = treatment),
    method = "wilcox.test",
    label = "p.format",
    size = 3
  ) +
  facet_wrap(~ response, ncol = 1) +
  scale_fill_manual(
    values = c("pre" = "#4D4D4D", "post" = "#1F78B4")
  ) +
  scale_color_manual(
    values = c("pre" = "#4D4D4D", "post" = "#1F78B4")
  ) +
  labs(
    x = NULL,
    y = "Fraction of CD4 T cells"
  ) +
  theme_nature
print(p2)

### 5.细胞成分Ro/e分析
library(Startrac)
library(circlize)
library(ComplexHeatmap)

CD4_pre <- subset(CD4, subset = treatment == "pre")
CD4_post <- subset(CD4, subset = treatment == "post")

Roe_pre <-calTissueDist(
  CD4_pre@meta.data,
  byPatient = F,
  colname.cluster = "celltype2", # 不同细胞分群，为结果中Ro/e矩阵的“行名”
  colname.patient = "sample_id", # 不同样本
  colname.tissue = "response", # 不同Group，为结果中Ro/e矩阵的“列名”
  method = "chisq", # 可选："chisq", "fisher", and "freq" 
  min.rowSum =0 # 整行的和为0则忽略
)
colnames(Roe_pre) <- paste0("pre_", colnames(Roe_pre))
Roe_post <-calTissueDist(
  CD4_post@meta.data,
  byPatient = F,
  colname.cluster = "celltype2", # 不同细胞分群，为结果中Ro/e矩阵的“行名”
  colname.patient = "sample_id", # 不同样本
  colname.tissue = "response", # 不同Group，为结果中Ro/e矩阵的“列名”
  method = "chisq", # 可选："chisq", "fisher", and "freq" 
  min.rowSum =0 # 整行的和为0则忽略
)
colnames(Roe_post) <- paste0("post_", colnames(Roe_post))

Roe <- cbind(Roe_pre, Roe_post)

breaks <- c(0.8, 0.9, 1.0, 1.1, 1.2)
colors <- c("#440154", "#31688e", "#35b779", "#fde725", "#ffffc0")
col_fun <- colorRamp2(breaks, colors)
# col_fun = colorRamp2(seq(0, 2.5, by = 0.5), c("#44045A", "#413E85", "#31698E", "#1F928B", "#35B777", "#91D540"))
pdf("analysis/4_functional_analysis/3_HSPC_Roe_heatmap.pdf", height = 5, width = 4)
Heatmap(as.matrix(Roe),
        show_heatmap_legend = T, 
        cluster_rows = T, 
        cluster_columns = T,
        row_names_side = 'right', 
        show_column_names = T,
        show_row_names = T,
        col = col_fun,
        row_names_gp = gpar(fontsize = 10),
        column_names_gp = gpar(fontsize = 10),
        heatmap_legend_param = list(
          title = "R_o/e",
          at =  seq(0.8, 1.2, by = 0.1),
          labels =  seq(0.8, 1.2, by = 0.1),
          legend_gp = gpar(fill = col_fun(seq(0.8, 1.2, by = 0.1)))
        ),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", Roe[i, j]), x, y, gp = gpar(fontsize = 8))
        }
      )
dev.off()

{
  library(scales)  # 用于颜色比例

# 1️⃣ 将矩阵转换为长格式
df <- as.data.frame(Roe)
df$Row <- rownames(df)
df_long <- df %>%
  pivot_longer(
    cols = -Row,
    names_to = "Column",
    values_to = "value"
  )

# 2️⃣ 自定义颜色渐变
breaks <- c(0.8, 0.9, 1.0, 1.1, 1.2)
colors <- c("#440154", "#31688e", "#35b779", "#fde725", "#ffffc0")

# 使用 scale_fill_gradientn 模拟 ComplexHeatmap 的 colorRamp2
col_fun <- scales::gradient_n_pal(colors, values = scales::rescale(breaks, to = c(0,1)))

# 3️⃣ 绘制热图
p <- ggplot(df_long, aes(x = Column, y = Row, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_gradientn(
    colors = colors,
    limits = c(min(breaks), max(breaks)),
    breaks = breaks,
    labels = breaks
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    panel.grid = element_blank()
  ) +
  labs(fill = "R_o/e")

# 4️⃣ 保存为 pdf
ggsave("analysis/4_functional_analysis/3_HSPC_Roe_heatmap_ggplot.pdf",
       p,
       width = 4, height = 5)
}


### 2.CFCT表达水平的比较 ----
set.seed(123)
## cell-level
{
  # 提取表达量数据
df <- FetchData(
  CD4,
  vars = c("CTCF", "IFNG", "celltype2", "sample_id", "response", "treatment")
)
df <- df %>%
  filter(celltype2 %in% paste0("C",1:15))
df$celltype2 <- factor(df$celltype2, levels = paste0("C",1:15))
df$treatment <- factor(df$treatment, levels = c("pre", "post"))
df$response  <- factor(df$response,  levels = c("NR", "R"))

# 循环计算所有组合的 CTCF 差异
celltypes <- paste0("C",1:15)
treatments <- c("pre", "post")

ctcf_stats <- expand.grid(
  ct = paste0("C",1:15),
  rs = c("NR", "R")
) %>%
  pmap_dfr(function(ct, rs) {

    cells_use <- WhichCells(
      CD4,
      expression = as.character(celltype2) == ct & as.character(response) == rs
    )

    sub <- subset(CD4, cells = cells_use)

    if (length(unique(sub$treatment)) < 2) return(NULL)

    de <- FindMarkers(
      sub,
      ident.1 = "post",
      ident.2 = "pre",
      group.by = "treatment",
      features = "CTCF",
      logfc.threshold = 0
    )

    tibble(
      celltype2 = ct,
      treatment = rs,
      avg_log2FC = de$avg_log2FC,
      p_val_adj = de$p_val_adj
    )
  })

ctcf_stats <- ctcf_stats %>%
  mutate(
    label = paste0(
      "log2FC = ", round(avg_log2FC, 2), "\n",
      "P_val_adj = ", signif(p_val_adj, 2)
    )
  )

theme_nature <- theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 0),
    axis.line = element_line(color = "black"),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "top",
    legend.title = element_blank()
  )

p_ctcf_cell <- ggplot(
  sample_frac(df, 0.3),
  aes(
    x = treatment,
    y = CTCF,
    fill = treatment
  )
) +
  geom_violin(
    outlier.shape = NA,
    alpha = 0.6,
    width = 0.65
  ) +
  geom_jitter(
    size = 0.1,
    alpha = 0.4,
    width = 0.15
  ) +
  facet_grid(
    response ~ celltype2,
    scales = "free_y"
  ) +
  # geom_text(
  #   data = ctcf_stats,
  #   aes(
  #     x = 1.5,
  #     y = Inf,
  #     label = label
  #   ),
  #   inherit.aes = FALSE,
  #   vjust = 1.2,
  #   size = 3
  # ) +
  stat_compare_means(
    comparisons = list(c("pre", "post")),
    method = "wilcox.test",
    label = "p.format",
    size = 3
  ) +
  scale_fill_manual(values = c("post" = "#B2182B", "pre" = "#2166AC")) +
  labs(
    x = NULL,
    y = "CTCF expression (log-normalized)"
  ) +
  theme_nature
print(p_ctcf_cell)
}
ggsave(filename = "results/3_analysis/01_CTCF_treatment.pdf", width = 10, height = 7)

## 不同癌种Treg中CTCF表达水平的比较
{
  treg_cells <- WhichCells(
  CD4,
  expression = celltype2 %in% c("C1", "C2", "C3", "C4", "C5")
)

Treg <- subset(CD4, cells = treg_cells)

df <- FetchData(
  Treg,
  vars = c("CTCF", "celltype2", "cancer_type", "sample_id", "response", "treatment")
)

sample_level_df <- df %>%
  dplyr::group_by(
    cancer_type,
    treatment,
    response,
    sample_id
  ) %>%
  dplyr::summarise(
    mean_CTCF = mean(CTCF, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

p <- ggplot(
  sample_level_df,
  aes(
    x = response,
    y = mean_CTCF,
    fill = response
  )
) +
  geom_violin(
    outlier.shape = NA,
    width = 0.65,
    alpha = 0.6
  ) +
  geom_jitter(
    aes(color = response),
    width = 0.15,
    size = 0.4,
    alpha = 0.4
  ) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.format",
    size = 3
  ) +
  facet_grid(
    treatment ~ cancer_type,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c("NR" = "#B2182B", "R" = "#2166AC")
  ) +
  scale_color_manual(
    values = c("NR" = "#B2182B", "R" = "#2166AC")
  ) +
  labs(
    x = NULL,
    y = "CTCF expression (log-normalized)",
    title = "CTCF expression in total Treg cells across cancer types"
  ) +
  theme_nature
print(p)
}

## CTCF相关性分析

clusters_use <- c("C1", "C2", "C3", "C4", "C5")
genes_use <- c("FOXP3", "IL2RA", "IKZF2", "TIGIT", "PDCD1", "CD40LG", "TregInst")
Idents(Treg) <- "celltype2"

# 计算 sample-level Spearman
spearman_by_cluster <- function(seu, cluster, genes, sample_col = "sample") {

  # 1. 选取 cluster 细胞
  cells <- WhichCells(seu, ident = cluster)

  # 2. 提取表达
  df <- FetchData(
    seu,
    vars = c("CTCF", genes, sample_col),
    cells = cells,
    layer = "data"
  )

  # 3. pseudobulk（sample-level mean）
  pb <- df %>%
    group_by(.data[[sample_col]]) %>%
    summarise(
      CTCF = mean(CTCF),
      across(all_of(genes), mean),
      .groups = "drop"
    )

  # 4. Spearman 相关
  res <- lapply(genes, function(g) {
    ct <- cor.test(pb$CTCF, pb[[g]], method = "spearman")
    data.frame(
      gene = g,
      rho = unname(ct$estimate),
      pval = ct$p.value
    )
  })

  res_df <- bind_rows(res)
  res_df$cluster <- cluster
  res_df
}

# 按Cluster批量计算
res_all <- bind_rows(
  lapply(clusters_use, function(clu) {
    spearman_by_cluster(
      seu = Treg,
      cluster = clu,
      genes = genes_use,
      sample_col = "sample_id"
    )
  })
)

# p值转星号
p_to_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}
res_all$gene <- factor(res_all$gene, levels = c("FOXP3", "IL2RA", "IKZF2", "TIGIT", "PDCD1", "CD40LG", "TregInst"))

res_all <- res_all %>%
  mutate(star = p_to_star(pval))

heatmap_plot <- ggplot(
  res_all,
  aes(x = cluster, y = gene, fill = rho)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(
    aes(label = star),
    color = "black",
    size = 5
  ) +
  scale_fill_gradient2(
    low = "#2166AC",   # 蓝：负相关
    mid = "white",
    high = "#B2182B",  # 红：正相关
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman rho"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.title = element_blank(),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )

heatmap_plot
ggsave(filename = "results/3_analysis/01_Treg_CTCF_correlation.pdf", height = 2.5, width = 4)

#### CTCF表达水平比较
# 1. 提取数据
ctcf_data <- FetchData(
  CD4,
  vars = c("CTCF", "dataset", "sample_id", "treatment", "response", "celltype2")
) %>%
  filter(celltype2 %in% c("C1", "C2", "C3", "C4", "C5", "C8", "C9", "C11", "C12", "C13"))

# 2. 过滤 + Sample-level 汇总
sample_level_data <- ctcf_data %>%
  filter(dataset == "GSE176021") %>%
  group_by(sample_id, treatment, response, celltype2) %>%
  summarise(
    CTCF = mean(CTCF, na.rm = TRUE),   # 也可以换成 median
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    treatment = factor(treatment, levels = c("pre", "post")),
    response = factor(response, levels = c("NR", "R"))
  )

ggplot(
  sample_level_data,
  aes(x = celltype2, y = CTCF, fill = response)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color = response),
    size = 0.5,
    alpha = 0.9,
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    )
  ) +
  stat_compare_means(
    aes(group = response),
    method = "wilcox.test",
    label = "p.format",
    size = 3
  ) +
  facet_wrap(~ treatment, ncol = 1) +
  scale_fill_manual(values = c("NR" = "#B2182B",
                               "R"  = "#2166AC")) +
  scale_color_manual(values = c("NR" = "#B2182B",
                                "R"  = "#2166AC")) +
  labs(
    x = "CD4 Subtype",
    y = "Mean CTCF Expression (sample-level)",
    title = "CTCF expression in CD4 cells across subtypes in GSE176021"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(filename = "results/3_analysis/02_GSE176021_CD4_CTCF_expression.pdf", height = 4, width = 12)

## Stability Score -- sample-level comparison
{
library(ggdist) # 核心库，推荐安装
library(ggsci)  # 高分期刊配色库
library(paletteer)
 
# 1. 提取数据
ctcf_data <- FetchData(
  CD4,
  vars = c("TregInst", "dataset", "sample_id", "treatment", "response", "celltype2")
) %>%
  filter(celltype2 %in% c("C1", "C2", "C3", "C4", "C5"))

# 2. 过滤 + Sample-level 汇总
sample_level_data <- ctcf_data %>%
  filter(dataset == "GSE176021") %>%
  group_by(sample_id, treatment, response, celltype2) %>%
  summarise(
    TregInst = mean(TregInst, na.rm = TRUE),   # 也可以换成 median
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    treatment = factor(treatment, levels = c("pre", "post")),
    response = factor(response, levels = c("NR", "R"))
  )

ggplot(sample_level_data, aes(x = celltype2, y = TregInst, fill = response)) +
  # 1. 使用半透明小提琴图背景
  geom_violin(aes(color = response), trim = FALSE, alpha = 0.1, 
              position = position_dodge(width = 0.8), scale = "width") +
  # 2. 箱线图变窄，去掉离群点，增加质感
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5,
               position = position_dodge(width = 0.8), color = "black") +
  # 3. 统计标注优化
  stat_compare_means(aes(group = response), method = "wilcox.test",
                     label = "p.signif", # CNS 常用星号标注，或者使用 "p.format"
                     label.y = max(sample_level_data$TregInst) * 1.1,
                     size = 3.5, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                    symbols = c("****", "***", "**", "*", "ns"))) +
  facet_wrap(~ treatment, ncol = 1) +
  # 4. 配色微调：Nature/Cell 常用更柔和但对比明确的颜色
  scale_fill_manual(values = c("NR" = "#E31A1C", "R" = "#1F78B4")) + 
  scale_color_manual(values = c("NR" = "#E31A1C", "R" = "#1F78B4")) +
  # 5. 主题细节提升
  theme_bw(base_size = 12) + 
  theme(
    panel.grid = element_blank(),           # 去掉网格线是 CNS 标配
    strip.background = element_blank(),     # 分面标签去框
    strip.text = element_text(face = "bold", size = 13),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",                # 图例置顶
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = expression(italic("TregInst") ~ "(sample-level)")) # 使用斜体标注基因或特定指标

ggsave(filename = "results/3_analysis/02_GSE176021_CD4_TregInst_score.pdf", height = 4, width = 4)
}


genesets <- list(
  PD1_perb = c(
    "FAM110A","MYL12A","TSPAN13","NRIP1","CIAPIN1","TRAF1","MT1M","DGAT1","ALDOA","KLRG1","ISG20",
    "PCYT1A","PLAC8","TRBC1","TRBC2","BHLHE40","NFKB2","SLC15A3","SYNGR2","CD74","GPR65","B2M",
    "CYBA","IKZF4","LGALS3BP","EBI3","IGKC","HOPX","CRBN","LIG1","VGLL4","IRF7","IER5","CCDC184",
    "DGAT2","YWHAQ","ODC1","RGS16","UAP1","DHRS3","TNFRSF18","TNFRSF9","TNFRSF4","ENO1"
  )
)
CD4 <- AddModuleScore_UCell(CD4, features = genesets, name = NULL, ncores = 12)

ctcf_data <- FetchData(
  CD4,
  vars = c("PD1_perb", "dataset", "sample_id", "treatment", "response", "celltype2")
) %>%
  filter(celltype2 %in% c("C1", "C2", "C3", "C4", "C5"))
sample_level_data <- ctcf_data %>%
  filter(dataset == "GSE176021") %>%
  group_by(sample_id, treatment, response, celltype2) %>%
  summarise(
    PD1_perb = mean(PD1_perb, na.rm = TRUE),   # 也可以换成 median
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    treatment = factor(treatment, levels = c("pre", "post")),
    response = factor(response, levels = c("NR", "R"))
  )
ggplot(sample_level_data, aes(x = celltype2, y = PD1_perb, fill = response)) +
  # 1. 使用半透明小提琴图背景
  geom_violin(aes(color = response), trim = FALSE, alpha = 0.1, 
              position = position_dodge(width = 0.8), scale = "width") +
  # 2. 箱线图变窄，去掉离群点，增加质感
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.5,
               position = position_dodge(width = 0.8), color = "black") +
  # 3. 统计标注优化
  stat_compare_means(aes(group = response), method = "wilcox.test",
                     label = "p.signif", # CNS 常用星号标注，或者使用 "p.format"
                     label.y = max(sample_level_data$PD1_perb) * 1.1,
                     size = 3.5, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                                                    symbols = c("****", "***", "**", "*", "ns"))) +
  facet_wrap(~ treatment, ncol = 1) +
  # 4. 配色微调：Nature/Cell 常用更柔和但对比明确的颜色
  scale_fill_manual(values = c("NR" = "#E31A1C", "R" = "#1F78B4")) + 
  scale_color_manual(values = c("NR" = "#E31A1C", "R" = "#1F78B4")) +
  # 5. 主题细节提升
  theme_bw(base_size = 12) + 
  theme(
    panel.grid = element_blank(),           # 去掉网格线是 CNS 标配
    strip.background = element_blank(),     # 分面标签去框
    strip.text = element_text(face = "bold", size = 13),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(color = "black"),
    legend.position = "top",                # 图例置顶
    legend.title = element_blank()
  ) +
  labs(x = NULL, y = expression(italic("PD1_perb_DN") ~ "(sample-level)")) # 使用斜体标注基因或特定指标
ggsave(filename = "results/3_analysis/02_GSE176021_CD4_PD1_perb_DN_score.pdf", height = 4, width = 4)


# 1. 提取数据
Ifng_data <- FetchData(
  CD4,
  vars = c("IFNG", "dataset", "sample_id", "treatment", "response", "celltype2")
) %>%
  filter(celltype2 %in% c("C1", "C2", "C3", "C4", "C5", "C8", "C9", "C11", "C12", "C13"))

# 2. 过滤 + Sample-level 汇总
sample_level_data <- Ifng_data %>%
  filter(dataset == "GSE176021") %>%
  group_by(sample_id, treatment, response, celltype2) %>%
  summarise(
    IFNG = mean(IFNG, na.rm = TRUE),   # 也可以换成 median
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    treatment = factor(treatment, levels = c("pre", "post")),
    response = factor(response, levels = c("NR", "R"))
  )

ggplot(
  sample_level_data,
  aes(x = celltype2, y = IFNG, fill = response)
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.6,
    alpha = 0.7,
    position = position_dodge(width = 0.75)
  ) +
  geom_jitter(
    aes(color = response),
    size = 0.5,
    alpha = 0.9,
    position = position_jitterdodge(
      jitter.width = 0.15,
      dodge.width  = 0.75
    )
  ) +
  stat_compare_means(
    aes(group = response),
    method = "wilcox.test",
    label = "p.format",
    size = 3
  ) +
  facet_wrap(~ treatment, ncol = 1) +
  scale_fill_manual(values = c("NR" = "#B2182B",
                               "R"  = "#2166AC")) +
  scale_color_manual(values = c("NR" = "#B2182B",
                                "R"  = "#2166AC")) +
  labs(
    x = "CD4 Subtype",
    y = "Mean IFNG Expression (sample-level)",
    title = "IFNG expression in CD4 cells across subtypes in GSE176021"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(filename = "results/3_analysis/02_GSE176021_CD4_IFNG_expression.pdf", height = 4, width = 12)

