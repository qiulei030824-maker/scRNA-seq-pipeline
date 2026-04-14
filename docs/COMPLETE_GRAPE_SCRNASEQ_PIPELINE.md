# 葡萄单细胞转录组完整分析流程
## 项目: PRJCA033302 (31个样本)

## 概述
本流程提供完整的葡萄(Vitis vinifera)单细胞转录组数据分析方案，针对10x Genomics平台数据，包含从原始FASTQ文件到质量控制过滤的完整流程。

## 文件清单

### 核心脚本
1. **`grape_scRNAseq_complete_pipeline.sh`** - 主流程脚本，整合所有步骤
2. **`seurat_grape_qc.R`** - Seurat质量控制R脚本
3. **`grape_scRNAseq_10x_analysis.py`** - Scanpy分析脚本（备用）

### 辅助脚本
4. **`process_all_grape_samples.sh`** - 批量处理脚本
5. **`README_grape_scRNAseq.md`** - 基础使用说明

## 完整分析流程

### 第一步：构建Cell Ranger参考基因组
**目标**: 创建高质量的葡萄核基因组参考

#### 1.1 下载基因组文件
```bash
# 从grapegenomics.com下载
# 需要下载的文件:
# 1. VvCabSauv_CC1.1.fasta - 基因组序列
# 2. VvCabSauv_CC1.1.gtf - 基因注释
# 保存到: /data5/qiulei/onepiece/data/reference/
```

#### 1.2 移除细胞器基因组
```bash
# 创建核基因组FASTA（只保留染色体序列）
grep -E "^>chr[0-9]+" VvCabSauv_CC1.1.fasta > VvCabSauv_CC1.1.nuclear.fasta

# 创建核基因组GTF（移除叶绿体和线粒体注释）
grep -v -E "chloroplast|mitochondria|ChrC|ChrM" VvCabSauv_CC1.1.gtf > VvCabSauv_CC1.1.nuclear.gtf

# 修复GTF中的gene_name（移除分号）
sed -i 's/;gene_name "/ gene_name "/g; s/";/"/g' VvCabSauv_CC1.1.nuclear.gtf
```

#### 1.3 构建Cell Ranger索引
```bash
cellranger mkref \
    --genome=VvCabSauv_CC1.1_nuclear \
    --fasta=VvCabSauv_CC1.1.nuclear.fasta \
    --genes=VvCabSauv_CC1.1.nuclear.gtf \
    --nthreads=16 \
    --memgb=64
```

### 第二步：运行Cell Ranger count
**目标**: 对每个样本进行基因表达定量

#### 2.1 单个样本处理
```bash
cellranger count \
    --id=CRR1417065 \
    --fastqs=/data5/qiulei/onepiece/data/raw/PRJCA033302 \
    --sample=CRR1417065 \
    --transcriptome=/data5/qiulei/onepiece/data/reference/VvCabSauv_CC1.1_nuclear \
    --include-introns \          # 必须添加，植物snRNA-seq需要
    --localcores=16 \
    --localmem=64 \
    --expect-cells=5000
```

#### 2.2 批量处理所有样本
使用提供的脚本：
```bash
# 修改脚本中的参数
chmod +x grape_scRNAseq_complete_pipeline.sh
./grape_scRNAseq_complete_pipeline.sh
```

### 第三步：Seurat质量控制与细胞过滤
**目标**: 过滤低质量细胞，准备下游分析

#### 3.1 运行Seurat分析
```bash
# 运行R脚本
Rscript seurat_grape_qc.R
```

#### 3.2 关键过滤参数
```r
# 葡萄数据专用过滤阈值
qc_thresholds <- list(
    nCount_RNA_min = 200,      # 最小UMI数
    nCount_RNA_max = 20000,    # 最大UMI数（防止doublets）
    nFeature_RNA_min = 180,    # 最小基因数
    percent_mt_max = 5,        # 最大线粒体基因百分比
    percent_chloro_max = 5     # 最大叶绿体基因百分比
)

# 葡萄基因模式
gene_patterns <- list(
    mito = "^MT-|^mt-|^MTRNR|^Mt-|^VIT_MT",      # 线粒体基因
    chloro = "^Cp-|^CHLORO|^CT-|^ChrC|^VIT_CP",  # 叶绿体基因
    ribo = "^RPS|^RPL|^Rps|^Rpl"                 # 核糖体基因
)
```

#### 3.3 输出文件
- `grape_filtered_seurat.rds` - 处理后的Seurat对象
- `grape_cell_metadata.csv` - 细胞元数据
- `grape_expression_matrix.mtx` - 基因表达矩阵
- `qc_violin_plot.png` - 质控小提琴图
- `qc_scatter_plots.png` - 质控散点图
- `grape_analysis_report.md` - 分析报告

## 环境要求

### 软件依赖
```bash
# Cell Ranger (10x Genomics)
# 下载: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest

# R包
install.packages(c("Seurat", "Matrix", "dplyr", "ggplot2", "patchwork"))
BiocManager::install(c("DoubletFinder", "scDblFinder", "harmony"))

# Python包 (可选)
pip install scanpy pandas numpy matplotlib scikit-learn
```

### 系统要求
- **内存**: 建议至少64GB RAM
- **存储**: 至少500GB可用空间
- **CPU**: 16核以上
- **操作系统**: Linux (推荐Ubuntu 20.04+)

## 执行步骤

### 快速开始
```bash
# 1. 下载基因组文件到指定目录
# 2. 运行完整流程
./grape_scRNAseq_complete_pipeline.sh

# 3. 检查输出
ls -la /data5/qiulei/onepiece/data/processed/VvPRJCA033302/
```

### 分步执行
```bash
# 步骤1: 构建参考基因组
./grape_scRNAseq_complete_pipeline.sh --step build-reference

# 步骤2: 运行Cell Ranger
./grape_scRNAseq_complete_pipeline.sh --step cellranger

# 步骤3: Seurat分析
./grape_scRNAseq_complete_pipeline.sh --step seurat
```

## 样本信息

### 31个样本列表
```
CRR1417065, CRR1417066, CRR1417067, CRR1417068, CRR1417069
CRR1417070, CRR1417071, CRR1417072, CRR1417073, CRR1417074
CRR1417075, CRR1417076, CRR1417077, CRR1417078, CRR1417079
CRR1417080, CRR1417081, CRR1417082, CRR1417083, CRR1417084
CRR1417085, CRR1417086, CRR1417087, CRR1417088, CRR1417089
CRR1417090, CRR1417091, CRR1417092, CRR1417093, CRR1417094
CRR1417095
```

### 数据目录结构
```
/data5/qiulei/onepiece/
├── data/
│   ├── raw/PRJCA033302/          # 原始FASTQ文件
│   │   ├── CRR1417065_r1.fastq.gz
│   │   ├── CRR1417065_r2.fastq.gz
│   │   └── ...
│   ├── reference/                 # 参考基因组
│   │   ├── VvCabSauv_CC1.1.fasta
│   │   ├── VvCabSauv_CC1.1.gtf
│   │   └── VvCabSauv_CC1.1_nuclear/  # Cell Ranger索引
│   └── processed/VvPRJCA033302/  # 处理结果
│       ├── CRR1417065/           # 每个样本的输出
│       ├── CRR1417066/
│       └── ...
└── logs/                         # 日志文件
```

## 质量控制要点

### 1. 植物数据特殊性
- **叶绿体基因**: 光合组织中含有大量叶绿体基因
- **线粒体基因**: 植物中线粒体基因前缀需要验证
- **核糖体基因**: 高表达可能指示低质量细胞

### 2. 过滤阈值调整建议
根据QC图调整：
- **nCount_RNA_max**: 如果发现明显的高UMI outlier，降低此值
- **percent_mt_max**: 如果大多数细胞<2%，可降低到3%
- **percent_chloro_max**: 根据组织类型调整，光合组织可放宽到10%

### 3. Doublet处理
- 使用scDblFinder预测doublets
- 典型doublet比例: 5-10%
- 移除doublets可提高聚类质量

## 故障排除

### 常见问题

#### Q1: Cell Ranger报错"Invalid GTF"
**原因**: GTF文件格式问题
**解决**: 检查并修复GTF文件中的gene_name列，移除分号

#### Q2: 内存不足
**解决**:
- 增加`--localmem`参数
- 分批处理样本
- 增加系统内存

#### Q3: Seurat加载数据失败
**解决**:
- 检查文件路径是否正确
- 确保Cell Ranger输出文件完整
- 检查R包版本兼容性

#### Q4: 线粒体基因识别错误
**解决**:
- 检查葡萄线粒体基因的实际前缀
- 修改`gene_patterns$mito`参数
- 查看GTF文件中的基因命名

### 日志检查
```bash
# 查看Cell Ranger日志
tail -f /data5/qiulei/onepiece/logs/cellranger_*.log

# 查看Seurat日志
Rscript seurat_grape_qc.R 2>&1 | tee seurat.log
```

## 下游分析建议

### 1. 降维与可视化
```r
# UMAP降维
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:30)

# t-SNE降维
filtered_seurat <- RunTSNE(filtered_seurat, dims = 1:30)

# 可视化
DimPlot(filtered_seurat, reduction = "umap", group.by = "sample")
```

### 2. 细胞聚类
```r
# 寻找邻居
filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:30)

# 细胞聚类
filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.5)

# 可视化聚类结果
DimPlot(filtered_seurat, reduction = "umap", label = TRUE)
```

### 3. 标记基因识别
```r
# 识别各簇的标记基因
markers <- FindAllMarkers(
    filtered_seurat,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
)

# 可视化标记基因
FeaturePlot(filtered_seurat, features = c("基因1", "基因2"))
```

### 4. 细胞类型注释
```r
# 基于已知标记基因注释
celltype_markers <- list(
    "表皮细胞" = c("基因A", "基因B"),
    "维管细胞" = c("基因C", "基因D"),
    "分生细胞" = c("基因E", "基因F")
)

# 细胞类型打分
filtered_seurat <- AddModuleScore(
    filtered_seurat,
    features = celltype_markers,
    name = "CellType"
)
```

## 性能优化

### 并行处理
```bash
# 使用GNU parallel并行处理样本
parallel -j 4 ./process_sample.sh {} ::: CRR1417065 CRR1417066 CRR1417067 CRR1417068
```

### 内存管理
- Cell Ranger: 设置`--localmem`根据可用内存调整
- Seurat: 使用稀疏矩阵节省内存
- 分批处理大样本

### 存储优化
- 使用压缩格式(.gz)存储中间文件
- 定期清理临时文件
- 使用符号链接管理大文件

## 验证与质控

### 质控检查清单
- [ ] Cell Ranger运行成功，输出文件完整
- [ ] 质控图显示合理的分布
- [ ] 过滤后细胞保留率>50%
- [ ] Doublet比例<10%
- [ ] 批次效应得到校正（多样本时）
- [ ] 高变基因数量合理(1000-3000)

### 数据质量指标
- **细胞数**: 每个样本>1000个高质量细胞
- **基因数**: 每个细胞>500个基因
- **线粒体百分比**: <5%
- **叶绿体百分比**: <5%（光合组织可放宽）
- **Doublet比例**: <10%

## 支持与联系

### 文档资源
- **Cell Ranger官方文档**: https://support.10xgenomics.com/
- **Seurat官方教程**: https://satijalab.org/seurat/
- **植物单细胞分析文献**: 见memory-bank中的参考文献

### 问题反馈
1. 检查日志文件中的错误信息
2. 验证输入文件格式和路径
3. 调整参数后重试
4. 查阅相关文档和文献

## 版本历史
- **v1.0 (2026-04-14)**: 初始版本，包含完整流程
- **v1.1 (2026-04-14)**: 添加故障排除和优化建议

## 许可证
本流程脚本遵循MIT许可证，可自由使用和修改。

---
**分析完成时间**: 2026-04-14  
**输出目录**: `/data5/qiulei/onepiece/data/processed/VvPRJCA033302/`  
**联系人**: Cline AI Assistant