library("Seurat")
library("ggplot2")
library(tidyverse)
library(gridExtra)
library(viridis)
library(Polychrome)
library(circlize)
library(ComplexHeatmap)

# Automatize it

lnames <- c("P1_Tumor","P10_Tumor","P2_Tumor","P3_Tumor","P4_Tumor","P5_Tumor","P6_Tumor",
            "P7_Tumor","P8_Tumor","P9_Tumor")

# Make the Seurat object with the raw data

dir <- getwd()
ldir <- list.files(paste0(dir,"/data_output/GSE144236/"), pattern = ".csv")

for(i in 1:length(ldir)){

p <- read.csv(paste0(dir,"./data_output/GSE144236/", ldir[i]), row.names = 1)

seurat <- CreateSeuratObject(counts = p, project = "221122", min.cells = 3, 
                             min.features = 10)
seurat@meta.data
setwd(dir)

# QC and filtering -------------------------------------------------------------

# calculate mitochondrial percentage

seurat$mitoPercent <- PercentageFeatureSet(seurat, pattern='^MT-')
head(seurat@meta.data) # we now have the mitopercent added to the metadata

# Quality check ----------------------------------------------------------------

# Just like Cell Ranger output, feature in the following results represents gene. 
# nFeature_RNA is the number of genes detected in each cell. nCount_RNA is the 
# total number of molecules detected within a cell. And each dot in the following 
# plots represents a cell.

p1 <- VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), 
              ncol = 3, pt.size = 0.0001)
ggsave(paste0(lnames[i],"_violin_plot_QC_unfiltered.png"), p1,"./figures/GSE144236/", width = 10,
       device = "png")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

p2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "mitoPercent") + 
  theme(legend.position="none")
p3 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
p2 + p3 # -0.44 and 0.90 show correlations between features and gene counts

# Here, we filter away cells that have unique feature counts(genes) over 5,000 
# or less than 200. We also filter away cells that have > 15% mitochondrial counts

seurat_f <- subset(seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & mitoPercent < 10)

p4 <- VlnPlot(seurat_f, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), 
              ncol = 3, pt.size = 0.0001)
print(p4) # adapt your threshold values to have good looking graphs
ggsave(paste0(lnames[i],"_violin_plot_QC_filtered.png"), p1,"./figures/GSE144236/", width = 10,
       device = "png")

p5 <- FeatureScatter(seurat_f, feature1 = "nCount_RNA", feature2 = "mitoPercent") + 
  theme(legend.position="none")
p6 <- FeatureScatter(seurat_f, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + 
  theme(legend.position="none")
p5 + p6 # -0.46 and 0.91 close to previous values, correlation 0.90 improved


## Data normalization-----------------------------------------------------------
# Normalized values are stored in seurat_f[["RNA"]]@data

seurat_f <- NormalizeData(seurat_f, normalization.method = "LogNormalize", 
                          scale.factor = 10000, verbose = FALSE)

rm(list = "seurat")
# Here we sample 10,000 reads counts from the large gene expression matrix to 
# visualize the gene expression distribution before and after normalization 
# separately (zeros are not included).

# set seed and put two plots in one figure
set.seed(1)
par(mfrow=c(1,2))
# original expression distribution
raw_geneExp = as.vector(seurat_f[['RNA']]@counts) %>% sample(10000)
raw_geneExp = raw_geneExp[raw_geneExp != 0]
hist(raw_geneExp)
# expression distribution after normalization
logNorm_geneExp = as.vector(seurat_f[['RNA']]@data) %>% sample(10000)
logNorm_geneExp = logNorm_geneExp[logNorm_geneExp != 0]
hist(logNorm_geneExp) # it is closer ta a normal distribution so 
# normalization was needed

# Normalize the entire data

seurat_f <- NormalizeData(seurat_f, normalization.method = "LogNormalize")

## Extract highly variable features --------------------------------------------


# We next calculate a subset of features that exhibit high cell-to-cell variation 
# in the dataset (i.e, they are highly expressed in some cells, and lowly 
#                 expressed in others). It is shown that (Brennecke et al. 2013) 
# focusing on these genes in downstream analysis helps to highlight biological 
# signal in single-cell datasets.
# 
# The procedure in Seurat is described in detail here (Stuart et al. 2019), and 
# improves on previous versions by directly modeling the mean-variance 
# relationship inherent in single-cell data, and is implemented in the 
# FindVariableFeatures() function. By default, Seurat returns 2,000 features 
# per dataset and these will be used in downstream analysis, like PCA.

seurat_f <- FindVariableFeatures(seurat_f, selection.method = "vst", nfeatures = 2000)
# need to check why vst is the chosen method

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seurat_f), 50)

# plot variable features with and without labels

par(mfrow=c(1,1))

p7 <- VariableFeaturePlot(seurat_f) + 
  theme(legend.position="top")
p8 <- LabelPoints(plot = p7, points = top10, repel = TRUE, max.overlaps = 50) + 
  theme(legend.position="top")
p8


## Scaling ---------------------------------------------------------------------

# Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing 
# step prior to dimensional reduction techniques like PCA. The ScaleData() function:
#   
# shifts the expression of each gene, so that the mean expression across cells is 0
# scales the expression of each gene, so that the variance across cells is 1
# this step gives equal weight in downstream analyses, so that highly-expressed 
# genes do not dominate the results of this are stored in seurat_f[["RNA"]]@scale.data

all_genes <- rownames(seurat_f)
seurat_f <- ScaleData(seurat_f, features = all_genes)

# Linear dimensional reduction--------------------------------------------------

seurat_f <- RunPCA(seurat_f, features = VariableFeatures(object = seurat_f))

DimPlot(seurat_f, reduction = "pca")

VizDimLoadings(seurat_f, dims = 1:2, reduction = "pca")

# In particular DimHeatmap() allows for easy exploration of the primary sources 
# of heterogeneity in a dataset, and can be useful when trying to decide which 
# PCs to include for further downstream analyses. Both cells and features are 
# ordered according to their PCA scores. Setting cells to a number plots the 
# ‘extreme’ cells on both ends of the spectrum, which dramatically speeds 
# plotting for large datasets. Though clearly a supervised analysis, we find 
# this to be a valuable tool for exploring correlated feature sets.

DimHeatmap(seurat_f, dims = 1, cells = 500, balanced = TRUE)

# To overcome the extensive technical noise in any single feature for scRNA-seq 
# data, Seurat clusters cells based on their PCA scores, with each PC essentially 
# representing a ‘metafeature’ that combines information across a correlated 
# feature set. The top principal components therefore represent a robust 
# compression of the dataset. Here, we choose first 20 PCs.

DimHeatmap(seurat_f, dims = 1:18, cells = 500, balanced = TRUE)
DimHeatmap(seurat_f, dims = 1:9, cells = 500, balanced = TRUE)

# Clustering the cells ---------------------------------------------------------

seurat_f <- FindNeighbors(seurat_f, dims = 1:14, verbose = FALSE)
seurat_f <- FindClusters(seurat_f, resolution = 0.5, verbose = FALSE)

# Look at cluster IDs of the first 5 cells
head(Idents(seurat_f), 5)


# UMAP / tsneq -----------------------------------------------------------------

seurat_f <- RunUMAP(seurat_f, dims = 1:20, verbose = FALSE)
DimPlot(seurat_f, reduction = "umap", label = T)

plot <- DimPlot(object = seurat_f)
LabelClusters(plot = plot, id = 'ident')
ggsave(paste0(lnames[i],"_UMAP_labelled.png"),last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png")

seurat_f <- RunTSNE(seurat_f, dims = 1:20, verbose = FALSE)
DimPlot(seurat_f, reduction = "tsne")
plot1 <- DimPlot(seurat_f, reduction = "tsne")
LabelClusters(plot1, id = "ident")

ggsave(paste0(lnames[i],"_tsne_labelled.png"),last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png")

# https://distill.pub/2016/misread-tsne/

saveRDS(seurat_f, file = paste0("./data_output/GSE144236/",lnames[i],"_seurat_f_processed.rds"))


## Exploiting gene markers------------------------------------------------------

# https://hbctraining.github.io/scRNA-seq/lessons/09_merged_SC_marker_identification.html

# Three objectives

## Determine gene markers for each cluster
## Identify cell types for each cluster
# Determine whether there’s a need to re-cluster based on cell type markers, 
# perhaps clusters need to be merged or split


# What are the cell type identities of clusters 7 and 20?
#   Do the clusters corresponding to the same cell types have biologically 
# meaningful differences? Are there subpopulations of these cell types?
#   Can we acquire higher confidence in these cell type identities by identifying 
# other marker genes for these clusters?
#   
#   There are a few different types of marker identification that we can explore 
# using Seurat to get to the answer of these questions. Each with their own 
# benefits and drawbacks:
#   
#   Identification of all markers for each cluster: this analysis compares each 
# cluster against all others and outputs the genes that are differentially 
# expressed/present.
# Useful for identifying unknown clusters and improving confidence in hypothesized 
# cell types.
# Identification of conserved markers for each cluster: This analysis looks for 
# genes that are differentially expressed/present within each condition first, 
# and then reports those genes that are conserved in the cluster across all 
# conditions. These genes can help to figure out the identity for the cluster.
# Useful with more than one condition to identify cell type markers that are 
# conserved across conditions.
# Marker identification between specific clusters: this analysis explores 
# differentially expressed genes between specific clusters.
# Useful for determining differences in gene expression between clusters that 
# appear to be representing the same celltype (i.e with markers that are similar) 
# from the above analyses

# Find markers for each clusters------------------------------------------------

# This type of analysis is typically recommended for when evaluating a single 
# sample group/condition. With the ` FindAllMarkers()` function we are comparing 
# each cluster against all other clusters to identify potential marker genes. 
# The cells in each cluster are treated as replicates, and essentially a 
# differential expression analysis is performed with some statistical test.
# By default, it's a wilcoxon rank sum test

markers <- FindAllMarkers(object = seurat_f,
                          only.pos = F,
                          logfc.threshold = 0.5)

write_csv(markers,paste0("./data_output/GSE144236/",lnames[i],"_markers_diff_per_cluster.csv"))

# Find conserved markers in each cluster (helps to identify celltype)-----------

# This function internally separates 
# out cells by sample group/condition, and then performs differential gene 
# expression testing for a single specified cluster against all other clusters 
# (or a second cluster, if specified). Gene-level p-values are computed for each 
# condition and then combined across groups using meta-analysis methods from the 
# MetaDE R package.
# 
# Before we start our marker identification we will explicitly set our default 
# assay, we want to use the original counts and not the integrated data.

# DefaultAssay(seurat_f) <- "RNA"
# 
# get_conserved <- function(cluster){
#   Seurat::FindConservedMarkers(seurat_f,
#                                ident.1 = cluster,
#                                only.pos = FALSE) %>%
#     rownames_to_column(var = "gene") %>%
#     left_join(y = unique(annotations[, c("gene_name", "description")]),
#               by = c("gene" = "gene_name")) %>%
#     cbind(cluster_id = cluster, .)
# }
# 
# FindConservedMarkers(seurat_f,
#                      ident.1 = "1",
#                      only.pos = FALSE)
# 
# 
# # Iterate function across desired clusters
# 
# conserved_markers <- map_dfr(c(0:18), get_conserved)

## Juliette and Chim-Kei's questions--------------------------------------------

# In which clusters are ROR2 etc expressed?

target <- c("ROR2","ROR1","STAT1","STAT3","CTNNB1")

FeaturePlot(object = seurat_f, 
            features = target,
            reduction = "tsne",
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE)

ggsave(paste0(lnames[i],"_ROR2_tSNE.png"),last_plot(),"./figures/GSE144236/", dpi = 500,
       width = 9, height = 9,
       device = "png")

VlnPlot(object = seurat_f, features = target)

ggsave(paste0(lnames[i],"_ROR2_violin_plot.png"),last_plot(),"./figures/GSE144236/", device = png, 
       dpi = 500, width = 9, height = 5)


# Markers for celltype identification

celltype <- c("KRT1","KRT5","KRT10","KRT14","LYZ","HLA-DRB1", "HLA-DRA","HLA-DQB2",
              "CD3D", "CD2","CD7","COL1A1", "COL1A2","LUM","MLANA", "DCT", "PMEL",
              "TFF3", "CLDN5", "VWF", "IGLL5", "IGJ"," MS4A1", "CD79A","CD163",
              "CD68","S100A8"," S100A9", "TREM1","CD1C", "CLEC10A","CLEC9A", "CADM1", 
              "XCR1", "CD207", "CD1A", "S100B","AXL", "IGFBP5", "PPP1R14A","CLEC4C",
              "IL3RA","CCR7", "CCL19","CD4","CD8")

dp <- DotPlot(seurat_f, features = celltype, cols = c("red","blue"))
dp

# Do it with complexheatmap

df <- dp$data

# the matrix for the scaled expression 
exp_mat<-df %>% 
  select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

row.names(exp_mat) <- exp_mat$features.plot  
exp_mat <- exp_mat[,-1] %>% as.matrix()

head(exp_mat)

## the matrix for the percentage of cells express a gene

percent_mat<-df %>% 
  select(-avg.exp, -avg.exp.scaled) %>%  
  pivot_wider(names_from = id, values_from = pct.exp) %>% 
  as.data.frame() 

row.names(percent_mat) <- percent_mat$features.plot  
percent_mat <- percent_mat[,-1] %>% as.matrix()

head(percent_mat)

Polychrome::swatch(viridis(20))

## any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])


cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}

reso <- 1200
length <- 3.25*reso/72

png(filename =paste0("./figures/GSE144236/heatmap",lnames[i],".png"))

plot <- Heatmap(exp_mat,
        heatmap_legend_param=list(title="expression"),
        column_title = "clustered dotplot", 
        col=col_fun,
        rect_gp = gpar(type = "none"),
        cell_fun = cell_fun,
        row_names_gp = gpar(fontsize = 5),
        row_km = 4,
        border = "black")
print(plot)
dev.off() 

}



