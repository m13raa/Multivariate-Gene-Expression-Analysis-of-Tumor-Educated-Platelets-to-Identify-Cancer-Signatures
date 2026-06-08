This project uses RNA sequencing data to explore whether gene expression patterns in blood platelets
can distinguish cancer from healthy samples and potentially differentiate between cancer types.

Tumor-educated platelets are platelets that have interacted with tumor cells and absorbed tumor-derived RNA,
acting as biomarkers for cancer detection. It serves as a useful took that allows for non-invasive testing
through blood samples instead of traditional biopsies. The key idea behind this project is that gene expression
profiles will reflect underlying disease states.

The dataset used in this project was downloaded from GEO and contains about 57,000 genes across 283 samples.
Out of these samples, 228 are identified as cancer, while 55 are identified as healthy samples.
The dataset includes multiple cancer types such as lung, breast, colorectal, etc, and makes for a high-dimensional
dataset with both continuous and categorical variables. 

The project focused on two main hypotheses: first, whether tumor-educated platelets can distinguish cancer from
healthy samples, and second, whether they can distinguish between different cancer types. To test this, a combination
of statistical and multivariate methods was used, such as DESeq2, PCA, ANOVA, clustering, and heatmaps.

Workflow:
The project began with filtering low-count genes to reduce noise. DESeq2 was used to normalize the data and perform
differential expression analysis. From the DEA, the significant genes were selected. PCA was performed to reduce
dimensionality and visualize variation within the data, and ANOVA testing was conducted to see whether the principal 
components differed between groups. Clustering was then applied to identify group structure, and heatmaps were used to
visualize gene expression patterns. 
