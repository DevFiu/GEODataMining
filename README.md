## **Differential Gene Expression Analysis**

**Project Introduction**

This project uses the limma package to analyze gene expression data and generate volcano plots and heatmaps.

**Data Source**

* GEO database: [https://www.ncbi.nlm.nih.gov/geo/](https://www.ncbi.nlm.nih.gov/geo/)

**Analysis Workflow**

1. Data preprocessing:
    * Read expression data and gene annotation files.
    * Filter missing values.
    * According to the gene annotation file, keep the probe row with the maximum value in each gene expression matrix.

2. Differential gene analysis:
    * Use the limma package to perform differential gene analysis.
    * Screen for significantly differentially expressed genes.

3. Volcano drawing:
    * Draw a volcano plot with logFC as the horizontal axis and -log10(P.Value) as the vertical axis.
    * Mark significantly differentially expressed genes.

4. Heatmap drawing:
    * Perform cluster analysis on significantly differentially expressed genes.
    * Draw a heatmap.

**Results**

* The differential gene analysis results are saved in the `./step2-deg.Rdata` file.
* The volcano plot is saved in the `./volcano.png` file.
* The heatmap is saved in the `./heatmap_top200_DEG.png` file.

**Running Environment**

* R language environment
* R packages: limma, dplyr, ggpubr, pheatmap

**Instructions**

1. Clone the project.
2. Install the required R packages.
3. Run the R script `./volcano.R`.

**Precautions**

* This project is for reference only, please adjust it according to the actual situation.
* The R language version and R package version may affect the results.

**Expected Results**

* Volcano:
    * Significantly up-regulated genes are located at the upper right of the volcano plot.
    * Significantly down-regulated genes are located at the lower left of the volcano plot.
* Heatmap:
    * The expression pattern of differentially expressed genes is clearly visible.
    * Genes with similar expression patterns are clustered together.

**Future Work**

* Further analysis of differentially expressed genes, such as functional annotation and pathway analysis.
