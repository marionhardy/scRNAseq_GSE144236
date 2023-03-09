# scRNAseq_GSE144236

scRNAseq analysis for SCBD lab based on GSE144236 dataset from the article:

Multimodal Analysis of Composition and Spatial Architecture in Human Squamous Cell Carcinoma.
Ji AL, Rubin AJ, Thrane K, Jiang S, Reynolds DL, Meyers RM, Guo MG, George BM, Mollbrink A, Bergenstråhle J, Larsson L, Bai Y, Zhu B, Bhaduri A, Meyers JM, Rovira-Clavé X, Hollmig ST, Aasi SZ, Nolan GP, Lundeberg J, Khavari PA. Cell. 2020 Sep 17;182(6):1661-1662. doi: 10.1016/j.cell.2020.08.043. PMID: 32946785. [ Free pmc article.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7391009/)

## Steps performed
### GSE file download

Single cell RNAseq from the 10 human patients (sCC + healthy) and the Cal27 and SCC13 xenografted cell lines

### Seurat object

- Creating the seurat object
- Quality check and filtering based on mitochondrial genes %, number of features and counts
- Normalization and scaling + regression of cell cycle genes (S and G2M phases), mitochondrial genes and counts
- PCA
- Clustering and visualization using both tSNE and UMAP
- Cluster identification based on DE markers using FindAllMarkers()
- Dotplot and clustering heatmap for cell type identification based on PBMC markers and author's markers
- Visualization of the features of interest 


### Data analysis










