**Differential Gene Expression Analysis**

**Introduction:**

This repository contains R code for performing differential gene expression analysis using data from the GEO database with accession number GSE233192.

**Data Preparation:**

The provided R script begins by clearing the workspace and loading necessary libraries such as limma and dplyr. It imports two datasets: "id_genesymbol.txt" containing gene IDs and symbols, and "GSE233192_series_matrix.txt" containing expression data. The gene symbols are cleaned to remove any empty values, and matching rows are selected from the expression dataset based on the gene IDs.

**Differential Expression Analysis:**

The script sets up experimental groups ('control' and 'treatment') and designs a linear model to analyze differential gene expression between these groups. Contrasts are created to compare the 'control' group with the 'treatment' group. Differential expression analysis is performed using the limma package, and significantly differentially expressed genes are identified and saved to a file named "step2-deg.Rdata".

**Visualization:**

Several plots are generated to visualize the results:
- A volcano plot displaying log-fold changes vs. significance levels.
- Scatter plots to visualize the relationship between log-fold changes and significance levels.
- A heatmap to visualize the expression profiles of the top 200 differentially expressed genes.

**Conclusion:**

This repository provides a comprehensive analysis pipeline for identifying and visualizing differentially expressed genes using R. Researchers can utilize this code to analyze gene expression data and gain insights into biological processes underlying experimental conditions.

**Data Source:**

The expression data used in this analysis can be accessed from the NCBI GEO database with accession number GSE233192.
