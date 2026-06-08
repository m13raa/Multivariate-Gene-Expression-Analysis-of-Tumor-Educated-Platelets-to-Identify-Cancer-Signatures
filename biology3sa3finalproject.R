#BIOLOGY 3SA3
#Multivariate Gene Expression Analysis of Tumor-Educated Platelets to Identify Cancer Signatures
#Workflow
#0. Load packages
#1. Bring in the data
#2. Sample group labels
#3. Use DESeq2 to a) normalize the data and b) perform differential expression analysis
#4. Keep only the significant genes (new matrix)
#5. Transpose the new matrix so we can run the PCA
#6. Run PCA
#7. Test PC1 and PC2 using lm() and anova()
#8. Run clustering on the PCA scores
#9. Make a heatmap

# Install DESeq2 through BiocManager if you don't have it already
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# Install the package for pheatmap
install.packages("pheatmap")

rm(list=ls()) #clear current environment

# SAIMIERAA's PART

# Load packages
library(tidyverse)
library(vegan) # allows for screeplot visualization
library(factoextra) # allows for pca visualization
library(DESeq2) # package not mentioned in lectures, allows for differential expression analysis and the normalization of the RNA-Seq data
library(pheatmap) # package not mentioned in lectures, allows for heatmap visualization

# Import data ( original csv has rows=genes, columns=samples)
ORIGINAL_TEP_DATA<-read.csv("biology_3sa3_final_project/biology_3sa3_final_project_data/GSE68086_TEP_data_matrix.csv", header=TRUE,check.names=FALSE) #I uploaded the files into my project folder so I can pull the files from there. Can change into file.choose() as well.
# Bring in the series file that identifies healthy vs. cancer samples and which type they are
ORIGINAL_SERIES_MATRIX<-read.csv("biology_3sa3_final_project/biology_3sa3_final_project_data/GSE68086_series_matrix.csv", header=TRUE,check.names=FALSE) #I uploaded the files into my project folder for easier access for this file as well

# Inspection
head(ORIGINAL_TEP_DATA)
str(ORIGINAL_TEP_DATA)
summary(ORIGINAL_TEP_DATA)

head(ORIGINAL_SERIES_MATRIX)
str(ORIGINAL_SERIES_MATRIX)
summary(ORIGINAL_SERIES_MATRIX)

# Creating a count matrix for format compatibility for DESeq2 and PCA
rownames(ORIGINAL_TEP_DATA) <- ORIGINAL_TEP_DATA[,1] # Turn the first column in the matrix (gene names) into row names
ORIGINAL_TEP_DATA <- ORIGINAL_TEP_DATA[,-1] # Take out the first column since its redundancy
ORIGINAL_TEP_DATA <- as.matrix(ORIGINAL_TEP_DATA) # Change into matrix so that the analyses can be run on the datafarame
mode(ORIGINAL_TEP_DATA) <- "numeric" # Ensure the data are comprised of integers for DESeq2
ORIGINAL_TEP_DATA<-round(ORIGINAL_TEP_DATA) # Round the data to the nearest integer for DESeq2
head(ORIGINAL_TEP_DATA) #Check first few rows of the count matrix

# Define the sample group labels and create the metadata frame to contain all information
snames<-colnames(ORIGINAL_TEP_DATA) #Get the sample names from the count matrix
sinfo<-data.frame(sname=snames) 
# Classify the samples based off the names given in the csv file using the grepl() function
sinfo$status_of_disease<- ifelse(grepl("HD|Control", sinfo$sname), "Healthy", "Cancer")
sinfo$type_of_cancer<- ifelse(grepl("HD|Control",sinfo$sname),"Healthy",
                              ifelse(grepl("Lung|NSCLC",sinfo$sname),"Lung Cancer",
                                     ifelse(grepl("Breast|BrCa",sinfo$sname),"Breast Cancer",
                                            ifelse(grepl("CRC",sinfo$sname),"Colorectal Cancer",
                                                   ifelse(grepl("GBM",sinfo$sname),"Glioblastoma",
                                                          ifelse(grepl("Panc|Pancr",sinfo$sname),"Pancreatic Cancer",
                                                                 ifelse(grepl("Liver|Chol",sinfo$sname),"Hepatobiliary", NA)))))))
# Inspection
head(sinfo)
summary(as.factor(sinfo$status_of_disease))
summary(as.factor(sinfo$type_of_cancer))

# Make the grouping variables turn into factors for DESeq2
sinfo$status_of_disease<-as.factor(sinfo$status_of_disease)
sinfo$type_of_cancer<-as.factor(sinfo$type_of_cancer)
rownames(sinfo)<-sinfo$sname #Have to match the metadata by making the sampe names equal to the row names

# Filter out any genes with reads less than 10 in 10 samples (at least)
keep_genes <- rowSums(ORIGINAL_TEP_DATA >= 10) >= 10
TEP_DATA_filtered <- ORIGINAL_TEP_DATA[keep_genes, ]

# Compare the dimensions of the unfiltered and filtered count matrix
dim(ORIGINAL_TEP_DATA)
dim(TEP_DATA_filtered)

# DEA with DESeq2 negative binomial test (accounts for over dispersion which is common in RNA Seq data)
# Compare expression of the cells based on the status_of_disease variable therefore Healthy vs. Cancer
dataset_DESeq2 <- DESeqDataSetFromMatrix(countData = TEP_DATA_filtered, colData = sinfo, design = ~ status_of_disease)
dataset_DESeq2 <- DESeq(dataset_DESeq2)
results_DESeq2 <- results(dataset_DESeq2, contrast = c("status_of_disease", "Cancer", "Healthy"))
dataframe_of_results<- as.data.frame(results_DESeq2)
dataframe_of_results$gene<- rownames(dataframe_of_results) #Add a column with the gene names 

# Check for significant genes (adjusted p-value < 0.01)
ordered_dataframe_of_results <- dataframe_of_results[order(dataframe_of_results$padj), ] #Lists the results by the order of adjusted p-value
head(dataframe_of_results)
summary(dataframe_of_results)

# We only want to keep the significant genes
sig_genes<- subset(ordered_dataframe_of_results, padj < 0.01 & abs(log2FoldChange)>1) #Subset the results to only keep significant genes with a log2 fold change greater than 1 or less than -1
# Inspection
dim(sig_genes)
head(sig_genes)
sig_gene_names <- sig_genes$gene
# Currently at 3809 genes, if we need to shorten further i made two seperate lists (one containing all, one containing top 200 for presentation)
top200_sig_gene_names <- sig_gene_names[1:min(200, length(sig_gene_names))] # Top 200 significant genes based on adjusted p-value and filtering through padj<0.01 and |log2|>1, selecting the first 200 genes
length(sig_gene_names)
length(top200_sig_gene_names)

# DEA summary table
DEAsummary<- data.frame(
  significant_genes = nrow(sig_genes),
  padj_cutoff = 0.01,
  abs_log2FC_cutoff = 1
)
DEAsummary

# Normalization of Data with DESeq2 results for PCA
stab_data <- vst(dataset_DESeq2, blind= FALSE) # Stabilized Variance
stab_matrix<- assay(stab_data) # Extract the normalized matrix from the DESeq2
stab_matrix_sig_genes <- stab_matrix[sig_gene_names, ] # Subset the normalized matrix to only include ALL the significant genes
stab_matrix_top200_sig_genes <- stab_matrix[top200_sig_gene_names, ] # Subset the normalized matrix to only include the top 200 significant genes

# Transpose the matrix to carry out PCA
PCAdata_all <- t(stab_matrix_sig_genes)
PCAdata_top200 <- t(stab_matrix_top200_sig_genes)

#JUSTINS PART 
#SAIMIERAA made edits to this section to ensure both hypotheses were being tested for All genes and Top 200 genes

#run PCA, with the data centered and scaled
pca_all <- prcomp(PCAdata_all, scale = TRUE)
pca_top200 <- prcomp(PCAdata_top200, scale = TRUE)

# Inspect PCA
summary(pca_all)
summary(pca_top200)       

# SCREE PLOTS - shows the variance explained by each principal component
screeplot(pca_all, bstick = TRUE, npcs = length(pca_all$sdev), type = "lines")
screeplot(pca_top200, bstick = TRUE, npcs = length(pca_top200$sdev), type = "lines")

# ANOVA - Test PC1 and PC2 using lm() and anova()
# Extract PC scores and have it merge with the metadata (match rownames)
# Hypothesis 1: Disease Status
pca_scores_all <- as.data.frame(pca_all$x)
pca_scores_all$sname <- rownames(pca_scores_all)
pca_scores_all <- merge(pca_scores_all, sinfo, by = "row.names")  # sinfo has rownames as sname                                                               
pca_scores_top200 <- as.data.frame(pca_top200$x)
pca_scores_top200$sname <- rownames(pca_scores_top200)
pca_scores_top200 <- merge(pca_scores_top200, sinfo, by = "row.names")

# All
lm_pc1_all <- lm(PC1 ~ status_of_disease, data = pca_scores_all)
anova(lm_pc1_all)
summary(lm_pc1_all)

lm_pc2_all <- lm(PC2 ~ status_of_disease, data = pca_scores_all)
anova(lm_pc2_all)
summary(lm_pc2_all)

# Top 200
lm_pc1_top200 <- lm(PC1 ~ status_of_disease, data = pca_scores_top200)
anova(lm_pc1_top200)
summary(lm_pc1_top200)

lm_pc2_top200 <- lm(PC2 ~ status_of_disease, data = pca_scores_top200)
anova(lm_pc2_top200)
summary(lm_pc2_top200)

# Hypothesis 2: Cancer Types - Remove the healthy samples for this analysis since they don't have a cancer type
pca_scores_allc <- subset(pca_scores_all, status_of_disease == "Cancer")
pca_scores_allc$type_of_cancer<- droplevels(pca_scores_allc$type_of_cancer) # Drop unused levels
pca_scores_top200c <- subset(pca_scores_top200, status_of_disease == "Cancer")
pca_scores_top200c$type_of_cancer<- droplevels(pca_scores_top200c$type_of_cancer) # Drop unused levels again

# All
lm_pc1_type_all <- lm(PC1 ~ type_of_cancer, data = pca_scores_allc)
anova(lm_pc1_type_all)
summary(lm_pc1_type_all)

lm_pc2_type_all <- lm(PC2 ~ type_of_cancer, data = pca_scores_allc)
anova(lm_pc2_type_all)
summary(lm_pc2_type_all)

# Top 200
lm_pc1_type_top200 <- lm(PC1 ~ type_of_cancer, data = pca_scores_top200c)
anova(lm_pc1_type_top200)
summary(lm_pc1_type_top200)

lm_pc2_type_top200 <- lm(PC2 ~ type_of_cancer, data = pca_scores_top200c)
anova(lm_pc2_type_top200)
summary(lm_pc2_type_top200)

# PCA SCATTER PLOTS

# simple (ALL genes) by disease status
plot(pca_scores_all$PC1, pca_scores_all$PC2,
     col = as.factor(pca_scores_all$status_of_disease),
     pch = 19,
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA Scatter Plot (All)")

legend("topright",
       legend = levels(pca_scores_all$status_of_disease),
       col = 1:length(levels(pca_scores_all$status_of_disease)),
       pch = 19) 
# simple (top 200) by disease status
plot(pca_scores_top200$PC1, pca_scores_top200$PC2,
     col = as.factor(pca_scores_top200$status_of_disease),
     pch = 19,
     xlab = "PC1",
     ylab = "PC2",
     main = "PCA Scatter Plot (Top 200)")

legend("topright",
       legend = levels(pca_scores_top200$status_of_disease),
       col = 1:length(levels(pca_scores_top200$status_of_disease)),
       pch = 19)

# biplot disease status (top 200)
fviz_pca_biplot(pca_top200,
                label = "none",   # removes ALL labels
                col.var = "black",  # arrows (genes)
                col.ind = pca_scores_top200$status_of_disease,
                palette = c("red", "blue"),
                addEllipses = FALSE,
                title = "PCA Scatter Plot by Disease Status (Top 200)")
# biplot cancer types (top 200)
fviz_pca_biplot(pca_top200,
                label = "none",   # removes ALL labels
                col.var = "black",  # arrows (genes)
                col.ind = pca_scores_top200$type_of_cancer,
                palette = c("red", "blue", "green", "purple", "orange","brown","cyan"),
                addEllipses = FALSE,
                title = "PCA Scatter Plot by Cancer Type (Top 200)")

# CLUSTERING - Select PCA scores for clustering
# We use principal components instead of raw gene data
# This reduces dimensionality and noise
# Here we take the first 5 PCs (can adjust if needed)

pc_data_cluster_top200 <- pca_scores_top200[, c("PC1", "PC2", "PC3", "PC4", "PC5")]

# Run hierarchical clustering using eclust()
# - Uses Ward’s method (minimizes variance within clusters)
# - Uses Euclidean distance (standard for PCA data)
# - k = 4 clusters (can adjust based on dendrogram)

hclust_top <- eclust(pc_data_cluster_top200,
                     FUNcluster = "hclust",
                     k = 4,
                     hc_method = "ward.D2",
                     hc_metric = "euclidean")

#  Plot the dendogram with cluster rectangles
# - rect = TRUE draws boxes around clusters
# - show_labels = FALSE removes clutter 
# - clean visualization of clustering structure

fviz_dend(hclust_top,
          rect = TRUE,
          show_labels = FALSE,
          cex = 0.6,
          main = "Hierarchical Clustering Dendrogram (Top 200)")

# Extract cluster assignments
# This gives each sample a cluster ID (1–4)

clusters_top <- hclust_top$cluster

# Add cluster labels back into PCA dataframe
pca_scores_top200$cluster <- clusters_top

# Evaluate clustering vs disease status
# This table shows how clusters correspond to:
# Healthy vs Cancer (main biological interpretation)

table(pca_scores_top200$cluster,
      pca_scores_top200$status_of_disease)

# Evaluate clustering vs cancer type
# More detailed breakdown of cluster composition

table(pca_scores_top200$cluster,
      pca_scores_top200$type_of_cancer)

# Colored dendrogram labels
# Highlights Healthy vs Cancer samples directly on dendrogram
#blue is healthy red is not

label_colors <- ifelse(pca_scores_top200$status_of_disease == "Healthy",
                       "blue", "red")

names(label_colors) <- rownames(pca_scores_top200)

fviz_dend(hclust_top,
          rect = TRUE,
          show_labels = TRUE,
          cex = 0.5,
          label_cols = label_colors,
          main = "Hierarchial Clustering Dendrogram Colored by Disease Status (Top 200)")

#try to color different cancer types and healthy (blue)

# Create factor (ensures consistent ordering)
cancer_factor <- as.factor(pca_scores_top200$type_of_cancer)

# Assign numeric colors
cancer_colors <- as.numeric(cancer_factor)

# Match labels to dendrogram
names(cancer_colors) <- rownames(pca_scores_top200)

fviz_dend(hclust_top,
          rect = TRUE,
          show_labels = TRUE,
          cex = 0.5,
          label_cols = cancer_colors,
          main = "Hierarchial Clustering Dendrogram Colored by Cancer Type (Top 200)")

legend("topright",
       inset = c(-0.15, 0),   # pushes legend outside plot
       legend = levels(cancer_factor),
       col = 1:length(levels(cancer_factor)),
       pch = 19,
       title = "Cancer Type")


# HEATMAP - visualize the expression of significant genes across samples
# Make a heatmap (of original sig genes expression, with sample clustering)
# Use pheatmap on normalized sig genes matrix (not transposed)
#all genes 
pheatmap(stab_matrix_sig_genes,
         scale = "row",  # Row z-score common for genes
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = sinfo[, c("status_of_disease", "type_of_cancer")],
         main = "Heatmap of Significant Genes (All)", )

# Ensure colnames(stab_matrix_top200_sig_genes) == rownames(sinfo)
annotation_col <- sinfo[colnames(stab_matrix_top200_sig_genes), c("status_of_disease", "type_of_cancer")]
print("Annotation dim:")
print(dim(annotation_col))  # Should match ncol(matrix) x 2

pheatmap(stab_matrix_top200_sig_genes,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Heatmap of Significant Genes (Top 200)",
         fontsize_row = 8)


top50_genes <- rownames(stab_matrix_top200_sig_genes)[1:50]
pheatmap(stab_matrix_top50_sig_genes <- stab_matrix_top200_sig_genes[top50_genes, ],
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = sinfo[colnames(stab_matrix_top50_sig_genes), c("status_of_disease", "type_of_cancer")],
         show_rownames = FALSE,
         show_colnames = FALSE,
         fontsize = 10,
         main = "Heatmap of Significant Genes (Top 50)")