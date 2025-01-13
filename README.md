# ATAC-Seq Analysis Scripts: Eosinophil Skin-Lung Crosstalk

This repository contains R scripts for analyzing ATAC-seq data as part of the study:

**"Bacterial skin infection primes bone marrow eosinophils to promote allergic skin-lung crosstalk."**

## Scripts Included

### 1. Data Preprocessing
- **`1_Peaks_Bam.R`**:  
  Processes BAM files, identifies consensus peaks, and creates count matrices.
  - **Output**: Peaks, annotated probes, and counts.

### 2. Differential Accessibility Analysis
- **`2_limma.R`**:  
  Runs differential accessibility analysis using `limma`.  
  - **Output**: Normalized data, volcano plots, significant hits.  
- **`2.1_annotated_volcano.R`**:  
  Creates annotated volcano plots of significant genes.  
- **`2.2_heatmap.R`**:  
  Generates heatmaps for top significant probes.  

### 3. Enrichment Analysis
- **`5_enricher.R`**:  
  Performs pathway enrichment analysis for chromatin states.
  - **Output**: Enriched pathways and visualizations.

## Requirements

- **R version**: 4.2.0 or higher  
- **R Libraries**:
  - `DiffBind`, `GenomicRanges`, `ChIPpeakAnno`, `limma`, `edgeR`
  - `tidyverse`, `ComplexHeatmap`, `pheatmap`, `enrichR`, `ggplot2`

## Usage

1. Upload your `meta.table.txt` and relevant input files to the same directory as the scripts.
2. Run scripts sequentially based on the pipeline described in the "Scripts Included" section.

## Contributors

- Mariem Radhouani
- Philipp Starkl  
- Natalia Nunes 
- Nikolaus Fortelny 

## Citation

Please cite the research:  
*"Bacterial skin infection primes bone marrow eosinophils to promote allergic skin-lung crosstalk."*
