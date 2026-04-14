#!/usr/bin/env Rscript
# ============================================================================
# 模块1: 数据加载与质量控制
# 功能: 加载单细胞数据，进行质控过滤，保存过滤后的数据
# 输入: Seurat对象文件 (.rds) 或 10X Genomics数据目录
# 输出: 过滤后的Seurat对象，质控图
# ============================================================================

# 加载必要的R包
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

# 设置随机种子
set.seed(42)

# 解析命令行参数
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("使用方法: Rscript 01_data_loading_qc.R <input_file> [output_dir]\n")
  cat("  <input_file>: Seurat对象文件(.rds)或10X数据目录\n")
  cat("  [output_dir]: 输出目录(可选，默认为当前目录)\n")
  quit(status = 1)
}

input_file <- args[1]
output_dir <- ifelse(length(args) > 1, args[2], ".")

# 创建输出目录
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
qc_plots_dir <- file.path(output_dir, "qc_plots")
dir.create(qc_plots_dir, recursive = TRUE, showWarnings = FALSE)

# 设置日志
log_file <- file.path(output_dir, "01_data_loading_qc.log")
sink(log_file, split = TRUE)

cat(paste(rep("=", 80), collapse = ""), "\n")
cat("模块1: 数据加载与质量控制\n")
cat(paste(rep("=", 80), collapse = ""), "\n\n")

# ============================================================================
# 步骤1: 数据加载
# ============================================================================

cat("步骤1: 数据加载\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

cat("输入文件:", input_file, "\n")

if (file.exists(input_file)) {
  if (grepl("\\.rds$", input_file, ignore.case = TRUE)) {
    # 加载Seurat对象
    cat("加载Seurat对象...\n")
    seurat_obj <- readRDS(input_file)
    cat("Seurat对象加载成功\n")
  } else if (dir.exists(input_file)) {
    # 加载10X Genomics数据
    cat("加载10X Genomics数据...\n")
    seurat_obj <- Read10X(data.dir = input_file)
    seurat_obj <- CreateSeuratObject(counts = seurat_obj, project = "scRNA")
    cat("10X数据加载成功\n")
  } else {
    stop("输入文件格式不支持。请提供.rds文件或10X数据目录")
  }
} else {
  stop("输入文件不存在: ", input_file)
}

# 检查数据
cat("\n数据基本信息:\n")
cat("细胞数:", ncol(seurat_obj), "\n")
cat("基因数:", nrow(seurat_obj), "\n")
cat("Assays:", names(seurat_obj@assays), "\n")

# 保存原始数据
original_file <- file.path(output_dir, "seurat_original.rds")
saveRDS(seurat_obj, original_file)
cat("原始数据已保存到:", original_file, "\n")

# ============================================================================
# 步骤2: 计算质控指标
# ============================================================================

cat("\n\n步骤2: 计算质控指标\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 计算线粒体基因百分比
cat("计算线粒体基因百分比...\n")

# 根据物种设置线粒体基因模式
species <- "arabidopsis"  # 默认拟南芥，可通过参数调整
if (species == "arabidopsis") {
  mt_pattern <- "^ATMG|^ATCG"
} else if (species == "human") {
  mt_pattern <- "^MT-"
} else if (species == "mouse") {
  mt_pattern <- "^mt-"
} else if (species == "grape") {
  mt_pattern <- "^VIT_MT|^MT-|^mt-"
} else {
  mt_pattern <- "^MT-|^mt-"
}

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = mt_pattern)
cat("线粒体基因百分比计算完成\n")

# 计算核糖体基因百分比
cat("计算核糖体基因百分比...\n")
ribo_pattern <- "^RPS|^RPL|^Rps|^Rpl"
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = ribo_pattern)
cat("核糖体基因百分比计算完成\n")

# 对于植物数据，计算叶绿体基因百分比
if (species %in% c("arabidopsis", "grape", "plant")) {
  cat("计算叶绿体基因百分比...\n")
  chloro_pattern <- "^ATCG|^VIT_CP|^Cp-|^CHLORO"
  seurat_obj[["percent.chloro"]] <- PercentageFeatureSet(seurat_obj, pattern = chloro_pattern)
  cat("叶绿体基因百分比计算完成\n")
}

# ============================================================================
# 步骤3: 质控可视化
# ============================================================================

cat("\n\n步骤3: 质控可视化\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 创建质控指标小提琴图
cat("生成质控小提琴图...\n")

qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
if (exists("percent.chloro", where = seurat_obj@meta.data)) {
  qc_features <- c(qc_features, "percent.chloro")
}
if (exists("percent.ribo", where = seurat_obj@meta.data)) {
  qc_features <- c(qc_features, "percent.ribo")
}

qc_violin <- VlnPlot(
  seurat_obj,
  features = qc_features,
  ncol = length(qc_features),
  pt.size = 0.1
)

# 保存小提琴图
violin_file <- file.path(qc_plots_dir, "qc_violin_plot.png")
ggsave(violin_file, qc_violin, width = 4 * length(qc_features), height = 8, dpi = 300)
cat("质控小提琴图已保存到:", violin_file, "\n")

# 创建散点图
cat("生成质控散点图...\n")

# nFeature_RNA vs nCount_RNA
scatter1 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  ggtitle("基因数 vs UMI数")

# percent.mt vs nCount_RNA
scatter2 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt") +
  ggtitle("线粒体百分比 vs UMI数")

# 如果有叶绿体数据，添加相应散点图
if (exists("percent.chloro", where = seurat_obj@meta.data)) {
  scatter3 <- FeatureScatter(seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.chloro") +
    ggtitle("叶绿体百分比 vs UMI数")
  combined_scatter <- scatter1 + scatter2 + scatter3
} else {
  combined_scatter <- scatter1 + scatter2
}

# 保存散点图
scatter_file <- file.path(qc_plots_dir, "qc_scatter_plots.png")
ggsave(scatter_file, combined_scatter, width = 12, height = 6, dpi = 300)
cat("质控散点图已保存到:", scatter_file, "\n")

# ============================================================================
# 步骤4: 细胞过滤
# ============================================================================

cat("\n\n步骤4: 细胞过滤\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 过滤前统计
cat("过滤前统计:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数:", nrow(seurat_obj), "\n")
cat("  平均UMI数:", mean(seurat_obj$nCount_RNA), "\n")
cat("  平均基因数:", mean(seurat_obj$nFeature_RNA), "\n")
cat("  平均线粒体百分比:", mean(seurat_obj$percent.mt), "\n")

if (exists("percent.chloro", where = seurat_obj@meta.data)) {
  cat("  平均叶绿体百分比:", mean(seurat_obj$percent.chloro), "\n")
}

# 设置过滤阈值
cat("\n应用过滤阈值:\n")

# 默认过滤阈值
filter_thresholds <- list(
  nFeature_RNA_min = 200,
  nFeature_RNA_max = 6000,
  nCount_RNA_min = 500,
  nCount_RNA_max = 30000,
  percent_mt_max = 10
)

# 根据物种调整阈值
if (species == "arabidopsis") {
  filter_thresholds$percent_mt_max <- 5
  filter_thresholds$percent_chloro_max <- 5
} else if (species == "grape") {
  filter_thresholds$percent_mt_max <- 5
  filter_thresholds$percent_chloro_max <- 10  # 葡萄光合组织可能含有较多叶绿体基因
}

cat("  最小基因数:", filter_thresholds$nFeature_RNA_min, "\n")
cat("  最大基因数:", filter_thresholds$nFeature_RNA_max, "\n")
cat("  最小UMI数:", filter_thresholds$nCount_RNA_min, "\n")
cat("  最大UMI数:", filter_thresholds$nCount_RNA_max, "\n")
cat("  最大线粒体百分比:", filter_thresholds$percent_mt_max, "\n")

if (exists("percent.chloro", where = seurat_obj@meta.data) && !is.null(filter_thresholds$percent_chloro_max)) {
  cat("  最大叶绿体百分比:", filter_thresholds$percent_chloro_max, "\n")
}

# 应用过滤
cat("\n应用细胞过滤...\n")

# 创建过滤条件
filter_condition <- (
  seurat_obj$nFeature_RNA > filter_thresholds$nFeature_RNA_min &
  seurat_obj$nFeature_RNA < filter_thresholds$nFeature_RNA_max &
  seurat_obj$nCount_RNA > filter_thresholds$nCount_RNA_min &
  seurat_obj$nCount_RNA < filter_thresholds$nCount_RNA_max &
  seurat_obj$percent.mt < filter_thresholds$percent_mt_max
)

if (exists("percent.chloro", where = seurat_obj@meta.data) && !is.null(filter_thresholds$percent_chloro_max)) {
  filter_condition <- filter_condition & (seurat_obj$percent.chloro < filter_thresholds$percent_chloro_max)
}

# 过滤细胞
filtered_seurat <- subset(seurat_obj, subset = filter_condition)

# 过滤后统计
cat("\n过滤后统计:\n")
cat("  保留细胞数:", ncol(filtered_seurat), "\n")
cat("  细胞保留率:", round(ncol(filtered_seurat) / ncol(seurat_obj) * 100, 2), "%\n")
cat("  保留基因数:", nrow(filtered_seurat), "\n")
cat("  平均UMI数:", mean(filtered_seurat$nCount_RNA), "\n")
cat("  平均基因数:", mean(filtered_seurat$nFeature_RNA), "\n")
cat("  平均线粒体百分比:", mean(filtered_seurat$percent.mt), "\n")

if (exists("percent.chloro", where = filtered_seurat@meta.data)) {
  cat("  平均叶绿体百分比:", mean(filtered_seurat$percent.chloro), "\n")
}

# ============================================================================
# 步骤5: 保存结果
# ============================================================================

cat("\n\n步骤5: 保存结果\n")
cat(paste(rep("-", 80), collapse = ""), "\n")

# 保存过滤后的Seurat对象
filtered_file <- file.path(output_dir, "seurat_filtered.rds")
saveRDS(filtered_seurat, filtered_file)
cat("过滤后的Seurat对象已保存到:", filtered_file, "\n")

# 保存细胞元数据
metadata_file <- file.path(output_dir, "cell_metadata.csv")
write.csv(filtered_seurat@meta.data, metadata_file)
cat("细胞元数据已保存到:", metadata_file, "\n")

# 保存过滤统计
stats_file <- file.path(output_dir, "qc_statistics.txt")
sink(stats_file)
cat("数据加载与质量控制统计\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")
cat("输入文件:", input_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("分析时间:", date(), "\n\n")

cat("过滤前统计:\n")
cat("  细胞数:", ncol(seurat_obj), "\n")
cat("  基因数:", nrow(seurat_obj), "\n")
cat("  平均UMI数:", round(mean(seurat_obj$nCount_RNA), 2), "\n")
cat("  平均基因数:", round(mean(seurat_obj$nFeature_RNA), 2), "\n")
cat("  平均线粒体百分比:", round(mean(seurat_obj$percent.mt), 2), "%\n")

if (exists("percent.chloro", where = seurat_obj@meta.data)) {
  cat("  平均叶绿体百分比:", round(mean(seurat_obj$percent.chloro), 2), "%\n")
}

cat("\n过滤后统计:\n")
cat("  保留细胞数:", ncol(filtered_seurat), "\n")
cat("  细胞保留率:", round(ncol(filtered_seurat) / ncol(seurat_obj) * 100, 2), "%\n")
cat("  保留基因数:", nrow(filtered_seurat), "\n")
cat("  平均UMI数:", round(mean(filtered_seurat$nCount_RNA), 2), "\n")
cat("  平均基因数:", round(mean(filtered_seurat$nFeature_RNA), 2), "\n")
cat("  平均线粒体百分比:", round(mean(filtered_seurat$percent.mt), 2), "%\n")

if (exists("percent.chloro", where = filtered_seurat@meta.data)) {
  cat("  平均叶绿体百分比:", round(mean(filtered_seurat$percent.chloro), 2), "%\n")
}

cat("\n过滤阈值:\n")
for (threshold_name in names(filter_thresholds)) {
  cat("  ", threshold_name, ":", filter_thresholds[[threshold_name]], "\n")
}
sink()

cat("质控统计已保存到:", stats_file, "\n")

# ============================================================================
# 完成
# ============================================================================

cat("\n\n" + paste(rep("=", 80), collapse = ""), "\n")
cat("模块1: 数据加载与质量控制 完成\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# 关闭日志
sink()

# 打印总结
cat("\n✅ 模块1执行完成!\n")
cat("📁 输出文件:\n")
cat("  - 过滤后的Seurat对象:", filtered_file, "\n")
cat("  - 细胞元数据:", metadata_file, "\n")
cat("  - 质控统计:", stats_file, "\n")
cat("  - 质控图目录:", qc_plots_dir, "\n")
cat("  - 日志文件:", log_file, "\n")

# 退出
quit(status = 0)