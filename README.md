# scRNAseq_GSE144236
## Data origin

scRNAseq analysis for SCBD lab based on GSE144236 dataset from the article:

Multimodal Analysis of Composition and Spatial Architecture in Human Squamous Cell Carcinoma.
Ji AL, Rubin AJ, Thrane K, Jiang S, Reynolds DL, Meyers RM, Guo MG, George BM, Mollbrink A, Bergenstråhle J, Larsson L, Bai Y, Zhu B, Bhaduri A, Meyers JM, Rovira-Clavé X, Hollmig ST, Aasi SZ, Nolan GP, Lundeberg J, Khavari PA. Cell. 2020 Sep 17;182(6):1661-1662. doi: 10.1016/j.cell.2020.08.043. [ PMID:32946785.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7391009/)

## Steps performed
### GSE file download

Single cell RNAseq from the 10 human patients (sCC + healthy) and the Cal27 and SCC13 xenografted cell lines

### Seurat object

- Creating the seurat object
- Quality check and filtering based on mitochondrial genes %, number of features and counts
- Normalization and scaling + regression of cell cycle genes (S and G2M phases), mitochondrial genes and counts
- PCA
- Clustering and visualization using both tSNE and UMAP
- Data integration to correct for batch effect
- Cluster identification based on DE markers using FindAllMarkers()
- Dotplot and clustering heatmap for cell type identification based on PBMC markers and author's markers
- Visualization of the features of interest 

### Data analysis

- Statistics about tumor composition
- Target gene expression level and identification of target gene positive and target gene negative celltypes and/or tumors
- Reclustering of cell types (particularly T cells and macrophages) to determine subpopulations profiles between groups of patients by determining the top10 varying genes per cell subtype
- Determining the top varying genes between epithelial ROR2+ and ROR2- patients 

### GSEA analysis

- Comparing ROR2+ and ROR2- epithelial cells, T cells and macrophages
- Comparing epithelial ROR2+ patients to epithelial ROR2- patients, checking epithelial , T cells and macrophages

### Gene expression correlation

- Using stat_cor() with pearson coefficient






