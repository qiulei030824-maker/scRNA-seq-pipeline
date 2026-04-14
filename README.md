# scRNA-seq Pipeline

Complete single-cell RNA-seq analysis pipeline with modular scripts for data loading, QC, normalization, clustering, annotation, and visualization.

## Overview

This repository contains a comprehensive single-cell RNA-seq analysis pipeline designed for plant and animal datasets. The pipeline includes:

- **Data preprocessing**: Quality control, filtering, and normalization
- **Dimensionality reduction**: PCA, UMAP, t-SNE
- **Clustering**: Cell type identification and marker gene discovery
- **Annotation**: Multiple annotation methods (AUCell, PCMaster, scPlantAnnotate)
- **Visualization**: Interactive plots and reports

## Pipeline Structure

### Core Analysis Modules

1. **01_data_loading_qc.R** - Data loading and quality control
2. **02_normalization_feature_selection.R** - Normalization and feature selection
3. **03_dimension_reduction_clustering.R** - Dimensionality reduction and clustering
4. **04_marker_gene_identification.R** - Marker gene identification
5. **05_aucell_annotation.R** - AUCell-based cell type annotation
6. **06_pcmaster_annotation.R** - PCMaster annotation method
7. **07_scplantannotate_annotation.R** - scPlantAnnotate for plant data
8. **08_integrated_scoring.R** - Integrated scoring system
9. **09_celltype_naming.R** - Final cell type naming

### Complete Workflows

- **grape_scRNAseq_complete_pipeline.sh** - Complete grape single-cell analysis pipeline
- **multi_species_scRNA_analysis_complete.R** - Multi-species analysis workflow
- **plant_scRNA_marker_selection_workflow.md** - Plant-specific marker selection guide

### Utility Scripts

- **download_fastq.sh** - Download FASTQ files from SRA
- **seurat_grape_qc.R** - Seurat-based quality control for grape data
- **scanpy_analysis.py** - Scanpy-based analysis pipeline
- **cell_annotation_pipeline.py** - Cell annotation pipeline

## Quick Start

### For Grape Data Analysis

```bash
# Run the complete grape analysis pipeline
chmod +x grape_scRNAseq_complete_pipeline.sh
./grape_scRNAseq_complete_pipeline.sh
```

### For General scRNA-seq Analysis

```bash
# Run the modular R pipeline
Rscript 01_data_loading_qc.R
Rscript 02_normalization_feature_selection.R
Rscript 03_dimension_reduction_clustering.R
# ... continue with other modules
```

## Installation

### Dependencies

#### R Packages
```r
install.packages(c("Seurat", "Matrix", "dplyr", "ggplot2", "patchwork", "harmony"))
BiocManager::install(c("AUCell", "SingleCellExperiment", "scater", "scran"))
```

#### Python Packages
```bash
pip install scanpy pandas numpy matplotlib scikit-learn leidenalg
```

#### System Tools
- Cell Ranger (for 10x Genomics data)
- STARsolo (alternative aligner)
- samtools, bedtools

## Usage Examples

### Example 1: Complete Grape Analysis
See [COMPLETE_GRAPE_SCRNASEQ_PIPELINE.md](COMPLETE_GRAPE_SCRNASEQ_PIPELINE.md) for detailed instructions.

### Example 2: Multi-species Analysis
```r
# Run the multi-species analysis
Rscript multi_species_scRNA_analysis_complete.R \
  --input data/expression_matrix.csv \
  --metadata data/cell_metadata.csv \
  --output results/
```

### Example 3: Cell Annotation
```python
# Run cell annotation pipeline
python cell_annotation_pipeline.py \
  --input filtered_cells.h5ad \
  --reference plant_cell_atlas.h5ad \
  --output annotated_cells.h5ad
```

## Data Requirements

### Input Formats
- **10x Genomics**: Cell Ranger output (filtered_feature_bc_matrix)
- **Expression matrices**: CSV, TSV, or H5AD format
- **Metadata**: CSV files with cell annotations

### Quality Control Metrics
- Minimum cells: > 1000 cells per sample
- Minimum genes: > 500 genes per cell
- Mitochondrial percentage: < 5%
- Doublet rate: < 10%

## Output Files

### Standard Outputs
- **QC reports**: HTML and PDF reports with quality metrics
- **Expression matrices**: Filtered and normalized matrices
- **Clustering results**: Cell cluster assignments
- **Marker genes**: Lists of cluster-specific markers
- **Annotation results**: Cell type annotations with confidence scores
- **Visualizations**: UMAP/t-SNE plots, heatmaps, violin plots

### Intermediate Files
- Raw and filtered Seurat objects (.rds)
- Scanpy AnnData objects (.h5ad)
- Log files and error reports

## Configuration

### Pipeline Parameters

Key parameters can be adjusted in configuration files:

1. **Quality control thresholds**: Modify in `01_data_loading_qc.R`
2. **Clustering resolution**: Adjust in `03_dimension_reduction_clustering.R`
3. **Annotation methods**: Configure in annotation scripts
4. **Visualization settings**: Customize in plotting functions

### Species-specific Settings

The pipeline includes species-specific configurations for:
- **Grape (Vitis vinifera)**: Special handling for chloroplast genes
- **Arabidopsis thaliana**: Plant-specific marker genes
- **Human/Mouse**: Standard mammalian settings

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce clustering resolution or subset data
2. **Annotation failures**: Check reference database compatibility
3. **Visualization issues**: Update ggplot2 or matplotlib versions
4. **Package conflicts**: Use virtual environments or containers

### Debugging Tips

- Check log files in the `logs/` directory
- Run scripts with `--verbose` flag for detailed output
- Test with small subsets before full analysis
- Validate input file formats

## Performance Optimization

### For Large Datasets
- Use sparse matrix representations
- Enable parallel processing where available
- Implement batch processing for multiple samples
- Use memory-efficient data structures

### Computational Resources
- **Memory**: 32GB+ RAM recommended
- **CPU**: 8+ cores for parallel processing
- **Storage**: 100GB+ for intermediate files
- **GPU**: Optional for accelerated computations

## Citation

If you use this pipeline in your research, please cite:

```
Cline AI Assistant. (2026). scRNA-seq Pipeline: Complete single-cell RNA-seq analysis workflow.
GitHub repository: https://github.com/qiulei030824-maker/scRNA-seq-pipeline
```

## License

This pipeline is released under the MIT License. See LICENSE file for details.

## Support

For questions, issues, or feature requests:
1. Check the documentation in the `docs/` directory
2. Review existing issues on GitHub
3. Contact the maintainers via GitHub issues

## Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with detailed description

## Acknowledgments

- 10x Genomics for Cell Ranger software
- Satija Lab for Seurat package
- Theis Lab for Scanpy package
- All contributors to the single-cell analysis ecosystem

---

**Last Updated**: 2026-04-14  
**Pipeline Version**: 1.0  
**Maintainer**: Cline AI Assistant