# GEODataMining

1. **Data Source**: The code begins with a comment indicating the source of the data, which is a GEO dataset with accession number GSE233192 from NCBI.

2. **Data Preparation**:
   - The code clears the workspace to start with a clean slate.
   - It loads the necessary R libraries: limma and dplyr.
   - It reads two files: "id_genesymbol.txt" containing gene IDs and symbols, and "GSE233192_series_matrix.txt" containing expression data.
   - It preprocesses the data by removing empty gene symbols and selecting corresponding expression data based on the gene IDs.

3. **Differential Gene Expression Analysis**:
   - The code sets up the experimental design, contrasts, and fits a linear model to the data using the limma package.
   - Differential expression analysis is performed using eBayes and topTable functions to identify significantly differentially expressed genes between two experimental conditions (control and treatment).
   - The results are saved to a file named "step2-deg.Rdata".

4. **Visualization**:
   - Several plots are generated to visualize the results:
     - A scatter plot of log-fold changes vs. -log10(p-values) is created to visualize differential expression, with genes color-coded based on significance and direction of change.
     - Another scatter plot is created to visualize the relationship between average expression and log-fold change, with genes color-coded based on significance levels.
     - A heatmap is generated to visualize the expression profiles of the top 200 differentially expressed genes.

Overall, the code performs data loading, preprocessing, differential expression analysis, and visualization using R libraries such as limma, dplyr, ggpubr, and pheatmap. It provides insights into gene expression changes between different experimental conditions.
