# Data source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE233192

# Clear workspace
rm(list = ls())

# Load required libraries
library(limma)  # For differential expression analysis
library(dplyr)  # For data manipulation

# Importing data
ids <- read.table("./id_genesymbol.txt", sep="\t", header=T) # Importing id_genesymbol table
exp <- read.table("./GSE233192_series_matrix.txt", sep="\t", header=T, row.names=1, comment.char="!") # Importing expression data table

# Data preprocessing
ids = ids[ids$GeneSymbol != '',]         # Removing empty values in the GeneSymbol column of ids
table(sort(table(ids$GeneSymbol)))       # Checking the distribution of GeneSymbol names in ids
length(unique(ids$GeneSymbol))           # Checking the number of unique GeneSymbol names in ids
exp = exp[rownames(exp) %in% ids$ID,]    # Selecting rows in exp with the same row names as IDs and saving to exp
ids = ids[match(rownames(exp), ids$ID),] # Selecting rows in ids with the same row names as IDs and saving to ids
head(ids)
head(exp)  # View data
length(unique(ids$GeneSymbol))
tmp = by(exp, ids$GeneSymbol, function(x) rownames(x) [which.max(rowMeans(x))]) # For each GeneSymbol in ids, keep the row with the maximum value among duplicates
probes = as.character(tmp)                            # Convert IDs to characters
exp = exp[rownames(exp) %in% probes,]                 # Keep data in exp with row names matching probes     
ids = ids[ids$ID %in% probes,]                        # Keep data in ids with row names matching probes
rownames(exp) = ids[match(rownames(exp), ids$ID), 4]  # Replace row names in exp with GeneSymbols from ids
save(exp, file = './step1-output.Rdata')

# Differential expression analysis
dat <- exp  
group_list = c(rep("control", times=4), rep("treat", times=4)) # Define groups based on website information: control and treat, each with four samples
design <- model.matrix(~0+factor(group_list))
colnames(design) = levels(factor(group_list)) 
contrast.matrix <- makeContrasts("control-treat", levels=design) # Contrast matrix: We want to compare control group with treat group for differential analysis

# Fit the model
fit <- lmFit(dat,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 

# Differential expression gene matrix
DiffEG <- topTable(fit2, coef=1, n=Inf) %>% na.omit()  
save(DiffEG, file='./step2-deg.Rdata')

# Plot Volcano
nrDEG = DiffEG
attach(nrDEG)
library(ggpubr)
df = nrDEG
df$v = -log10(P.Value) # Add a column 'v' to df with values as -log10(P.Value)
ggscatter(df, x = "logFC", y = "v", size=0.5) # Scatter plot
df$col = ifelse(df$P.Value>0.01, 'stable',    # If the P.Value for a gene is >0.01, it's labeled as 'stable'.
                ifelse( df$logFC >2, 'up',    # Else, if P.Value<0.01, further check logFC: if >2, labeled as 'up' (up-regulated),
                        ifelse( df$logFC < -2, 'down', 'stable'))) # if < -2, labeled as 'down' (down-regulated), else 'stable'.
table(df$col)
df$name = rownames(df)
ggscatter(df, x="logFC", y="v", size=0.5, color='col')
ggscatter(df, x="logFC", y="v", color="col", size=0.5, label="name", repel=T,
   #label.select = rownames(df)[df$g != 'stable'] ,
   label.select = head(rownames(df)),  # Select some genes to display on the plot
   palette = c("#00AFBB", "#E7B800", "#FC4E07") )
ggsave('./volcano.png', width=8, height=6) # Save plot

ggscatter(df, x="AveExpr", y="logFC", size=0.2)
df$p_c = ifelse(df$P.Value<0.001, 'p<0.001', ifelse(df$P.Value<0.01, '0.001<p<0.01', 'p>0.01'))
table(df$p_c)
ggscatter(df, x="AveExpr", y="logFC", color="p_c", size=0.2, palette=c("green", "red", "black") )
ggsave(df, './MA.png', width=8, height=6)

# Plot Heatmap
library(pheatmap)
x = DiffEG$logFC             # Extract the logFC column from DiffEG and assign it to x
names(x) = rownames(DiffEG)  # Assign row names of DiffEG as names to 
cg = c(names(head(sort(x), 100)), names(tail(sort(x), 100))) # Sort x in ascending order, take the top 100 and bottom 100, and assign their corresponding probe names to cg
pheatmap(dat[cg,], show_colnames=F, show_rownames=F) # Draw heatmap using the matrix obtained by selecting rows of dat based on cg
n = t(scale(t(dat[cg,]))) # Normalize the log-ratio values using "scale", as it is applied across different groups, transpose the dat to match the required format
                          # since the scale function is applied across different groups, it requires row names to be samples, hence the need to transpose with t(dat[cg,])
n[n>2] = 2
n[n< -2] = -2
n[1:4,1:4]
pheatmap(n,show_colnames=F, show_rownames=F) 
ac = data.frame(g=group_list)
rownames(ac) = colnames(n)  # Assign the row names of ac, which represent group information, to the column names of n, representing group information in the heatmap
pheatmap(n, show_colnames=F, show_rownames=F, cluster_cols=F, annotation_col=ac, filename='./heatmap_top200_DEG.png')  # Column annotation information is ac, i.e., group information
