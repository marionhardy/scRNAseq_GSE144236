
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(viridis)
library(Polychrome)
library(circlize)
library(ComplexHeatmap)
library(patchwork)

## Get subset of patients expressing ROR2 or not (P5, P8)

sptROR2 = subset(seurat_all, subset = orig.ident %in% c("P1","P2","P4","P6","P7",
                                                        "P9","P10")&condition=="Tumor")

sptnoROR2 = subset(seurat_all, subset = orig.ident %in% c("P3","P5","P8")&condition=="Tumor")

# Re cluster myeloid populations and CD4 CD8+ cells
# ROR2neg epithelial cells patients

sptnoROR2_mt = subset(sptnoROR2, subset = paper_clusters %in% c("Myeloid cells","T cells"))

sptnoROR2_mt <- FindNeighbors(sptnoROR2_mt, dims = 1:14)
sptnoROR2_mt <- FindClusters(sptnoROR2_mt, resolution = 2, verbose = FALSE)

sptnoROR2_mt <- RunUMAP(sptnoROR2_mt, dims = 1:20, verbose = FALSE)
DimPlot(sptnoROR2_mt, reduction = "umap", label = T)

plot <- DimPlot(object = sptnoROR2_mt)
LabelClusters(plot = plot, id = 'ident')

ggsave("myeloid_tcells_patient_tumor_ror2neg_UMAP_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png", height = 7, width = 9)


sptnoROR2_mt <- RunTSNE(sptnoROR2_mt, dims = 1:20, verbose = FALSE)
DimPlot(sptnoROR2_mt, reduction = "tsne")
plot1 <- DimPlot(sptnoROR2_mt, reduction = "tsne")
LabelClusters(plot1, id = "ident")
ggsave("myeloid_tcells_patient_tumor_ror2neg_res2_tSNE_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png", height = 7, width = 9)


Idents(sptnoROR2_mt) <- sptnoROR2_mt@meta.data$paper_clusters
plot1 <- DimPlot(sptnoROR2_mt, reduction = "tsne")
LabelClusters(plot1, id = "ident")
ggsave("myeloid_tcells_patient_tumor_ror2neg_tSNE_author_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png", height = 7, width = 9)

Idents(sptnoROR2_mt) <- sptnoROR2_mt@meta.data$attr_clusters
plot1 <- DimPlot(sptnoROR2_mt, reduction = "tsne")
LabelClusters(plot1, id = "ident")
ggsave("myeloid_tcells_patient_tumor_ror2neg_tSNE_attributed_labelled.png",last_plot(),
       "./figures/GSE144236/", dpi = 500 ,
       device = "png", height = 7, width = 9)

# ROR2pos epithelial patients

sptROR2_mt = subset(sptROR2, subset = paper_clusters %in% c("Myeloid cells","T cells"))

sptROR2_mt <- FindNeighbors(sptROR2_mt, dims = 1:14)
sptROR2_mt <- FindClusters(sptROR2_mt, resolution = 2, verbose = FALSE)

sptROR2_mt <- RunUMAP(sptROR2_mt, dims = 1:20, verbose = FALSE)
DimPlot(sptROR2_mt, reduction = "umap", label = T)

plot <- DimPlot(object = sptROR2_mt)
LabelClusters(plot = plot, id = 'ident')

ggsave("myeloid_tcells_patient_tumor_ror2pos_UMAP_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png", height = 7, width = 9)


sptROR2_mt <- RunTSNE(sptROR2_mt, dims = 1:20, verbose = FALSE)
DimPlot(sptROR2_mt, reduction = "tsne")
plot1 <- DimPlot(sptROR2_mt, reduction = "tsne")
LabelClusters(plot1, id = "ident")
ggsave("myeloid_tcells_patient_tumor_ror2pos_res2_tSNE_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png", height = 7, width = 9)


Idents(sptROR2_mt) <- sptROR2_mt@meta.data$paper_clusters
plot1 <- DimPlot(sptROR2_mt, reduction = "tsne")
LabelClusters(plot1, id = "ident")
ggsave("myeloid_tcells_patient_tumor_ror2pos_tSNE_author_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png", height = 7, width = 9)

Idents(sptROR2_mt) <- sptROR2_mt@meta.data$attr_clusters
plot1 <- DimPlot(sptROR2_mt, reduction = "tsne")
LabelClusters(plot1, id = "ident")
ggsave("myeloid_tcells_patient_tumor_ror2pos_tSNE_attributed_labelled.png",last_plot(),
       "./figures/GSE144236/", dpi = 500 ,
       device = "png", height = 7, width = 9)

# Markers for celltype identification

celltype <- c("KRT1","KRT5","KRT10","KRT14","LYZ","HLA-DRB1", "HLA-DRA","HLA-DQB2",
              "CD3D", "CD2","CD7","COL1A1", "COL1A2","LUM","MLANA", "DCT", "PMEL",
              "TFF3", "CLDN5", "VWF", "IGLL5", "IGJ"," MS4A1", "CD79A","CD163",
              "CD68","S100A8"," S100A9", "TREM1","CD1C", "CLEC10A","CLEC9A", "CADM1", 
              "XCR1", "CD207", "CD1A", "S100B","AXL", "IGFBP5", "PPP1R14A","CLEC4C",
              "IL3RA","CCR7", "CCL19","CD4","FOXP3","CD8A","MKI67","XCL1","CD44","LY6C",
              "CCT6A","CD45RA","CD45RO")

Idents(sptROR2_mt) <- sptROR2_mt@meta.data$RNA_snn_res.2

dp <- DotPlot(sptROR2_mt, features = celltype, cols = c("red","blue"))
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

png(filename ="./figures/GSE144236/myeloid_tcells_res2_ROR2p_newclusters_heatmap.png")

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




## Recluster the macrophages----------------------------------------------------

macr <- c("CD80","CD86","CCL22","CD64","CD32","CD40","SOCS3","SOCS1","IDO",
              "CD206","CD163","MRC1","CD68","CD23","TGM2","VEGF")

# sptROR2 = subset(seurat_all, subset = orig.ident %in% c("P1","P2","P4","P6","P7",
#                                                         "P9","P10")&condition=="Tumor")

sptROR2_macr = subset(sptROR2, subset = attr_clusters %in% c("Macrophage_1",
                                                             "Macrophages_2"))

sptROR2_macr <- FindNeighbors(sptROR2_macr, dims = 1:14)
sptROR2_macr <- FindClusters(sptROR2_macr, resolution = 0.5, verbose = FALSE)

sptROR2_macr <- RunTSNE(sptROR2_macr, dims = 1:14, verbose = FALSE)
DimPlot(sptROR2_macr, reduction = "tsne")
plot1 <- DimPlot(sptROR2_macr, reduction = "tsne")
LabelClusters(plot1, id = "ident")
ggsave("macr_patient_tumor_ror2pos_res05_tSNE_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 300 ,
       device = "png", height = 5, width = 7)

Idents(sptROR2_macr) <- sptROR2_macr@meta.data$RNA_snn_res.0.5

dp <- DotPlot(sptROR2_macr, features = macr, cols = c("red","blue"))
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

## any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])


cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}

png(filename ="./figures/GSE144236/macr_res2_ROR2p_newclusters_heatmap.png")

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

## Recluster the T cells--------------------------------------------------------

tcell = c("CCR7","FOXP3","CD4", "CD8A","CD44","LY6C","CCT6A","CD45RO","CD45RA")

sptROR2_tcells = subset(sptROR2, subset = attr_clusters == "T cells")

sptROR2_tcells <- FindNeighbors(sptROR2_tcells, dims = 1:14)
sptROR2_tcells <- FindClusters(sptROR2_tcells, resolution = 0.5, verbose = FALSE)

sptROR2_tcells <- RunTSNE(sptROR2_tcells, dims = 1:14, verbose = FALSE)
DimPlot(sptROR2_tcells, reduction = "tsne")
plot1 <- DimPlot(sptROR2_tcells, reduction = "tsne")
LabelClusters(plot1, id = "ident")
ggsave("tcells_patient_tumor_ror2pos_res05_tSNE_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 300 ,
       device = "png", height = 5, width = 7)


Idents(sptROR2_tcells) <- sptROR2_tcells@meta.data$RNA_snn_res.0.5

dp <- DotPlot(sptROR2_tcells, features = tcell, cols = c("red","blue"))
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

## any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])


cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}

png(filename ="./figures/GSE144236/tcells_res2_ROR2p_newclusters_heatmap.png")

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


### Now on patients with no epithelial ROR2-------------------------------------

## Recluster the macrophages----------------------------------------------------

macr <- c("CD80","CD86","CCL22","CD64","CD32","CD40","SOCS3","SOCS1","IDO",
          "CD206","CD163","MRC1","CD68","CD23","TGM2","VEGF","XCL1")

# sptnoROR2 = subset(seurat_all, subset = orig.ident %in% c("P5","P8")&condition=="Tumor")

table(sptnoROR2@meta.data$attr_clusters)

sptnoROR2_macr = subset(sptnoROR2, subset = attr_clusters %in% c("Macrophage_1",
                                                             "Macrophages_2"))

Idents(sptnoROR2)= "orig.ident"
sptnoROR2_macr <- FindNeighbors(sptnoROR2_macr, dims = 1:14)
sptnoROR2_macr <- FindClusters(sptnoROR2_macr, resolution = 0.5, verbose = FALSE)

sptnoROR2_macr <- RunTSNE(sptnoROR2_macr, dims = 1:14, verbose = FALSE)
DimPlot(sptnoROR2_macr, reduction = "tsne")
plot1 <- DimPlot(sptnoROR2_macr, reduction = "tsne")
LabelClusters(plot1, id = "ident")
ggsave("macr_patient_tumor_ror2neg_res05_tSNE_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 300 ,
       device = "png", height = 5, width = 7)

Idents(sptnoROR2_macr) <- sptnoROR2_macr@meta.data$RNA_snn_res.0.5

dp <- DotPlot(sptnoROR2_macr, features = macr, cols = c("red","blue"))
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

## any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])


cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}

png(filename ="./figures/GSE144236/macr_res2_ROR2n_newclusters_heatmap.png")

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

## Recluster the T cells--------------------------------------------------------

tcell = c("CCR7","FOXP3","CD4", "CD8A","CD44","LY6C","CCT6A","CD45RO","CD45RA", "XCL1")

sptnoROR2_tcells = subset(sptnoROR2, subset = attr_clusters == "T cells")

sptnoROR2_tcells <- FindNeighbors(sptnoROR2_tcells, dims = 1:14)
sptnoROR2_tcells <- FindClusters(sptnoROR2_tcells, resolution = 0.5, verbose = FALSE)

sptnoROR2_tcells <- RunTSNE(sptnoROR2_tcells, dims = 1:14, verbose = FALSE)
DimPlot(sptnoROR2_tcells, reduction = "tsne")
plot1 <- DimPlot(sptnoROR2_tcells, reduction = "tsne")
LabelClusters(plot1, id = "ident")
ggsave("tcells_patient_tumor_ror2neg_res05_tSNE_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 300 ,
       device = "png", height = 5, width = 7)


Idents(sptnoROR2_tcells) <- sptnoROR2_tcells@meta.data$RNA_snn_res.0.5

dp <- DotPlot(sptnoROR2_tcells, features = tcell, cols = c("red","blue"))
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

## any value that is greater than 2 will be mapped to yellow
col_fun = circlize::colorRamp2(c(-1, 0, 2), viridis(20)[c(1,10, 20)])


cell_fun = function(j, i, x, y, w, h, fill){
  grid.rect(x = x, y = y, width = w, height = h, 
            gp = gpar(col = NA, fill = NA))
  grid.circle(x=x,y=y,r= percent_mat[i, j]/100 * min(unit.c(w, h)),
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))}

png(filename ="./figures/GSE144236/tcell_res2_ROR2n_newclusters_heatmap.png")

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



## Assigning identities to all the clusters

sptROR2_macr@meta.data$sub_clust05 = NA

sptROR2_macr@meta.data <-
  sptROR2_macr@meta.data %>%
  mutate(sub_clust05 =
           ifelse(seurat_clusters == "0", "M1 macrophages_1",
                         ifelse(seurat_clusters == "1", "M0 macrophages",
                                ifelse(seurat_clusters == "2", "M1 macrophages_2",
                                       ifelse(seurat_clusters == "3", "M1/M2 macrophages",
                                              ifelse(seurat_clusters == "4", "M2 macrophages", 
                                                            NA))))))

Idents(sptROR2_macr) <- sptROR2_macr@meta.data$sub_clust05
temp <- DimPlot(sptROR2_macr, reduction = "tsne")
LabelClusters(temp, id = "ident")

ggsave("macr_patient_tumor_ror2p_tsne_subclust.png",last_plot(),
       "./figures/GSE144236/", device = png, dpi = 300,
       width = 8, height = 5)

saveRDS(spt, file = "./data_output/GSE144236/norm_log_scale/spt.rds")
saveRDS(sptROR2, file = "./data_output/GSE144236/norm_log_scale/sptROR2.rds")
saveRDS(sptnoROR2, file = "./data_output/GSE144236/norm_log_scale/sptnoROR2.rds")
saveRDS(sptROR2_macr, file = "./data_output/GSE144236/norm_log_scale/sptROR2_macr.rds")
saveRDS(sptROR2_tcells, file = "./data_output/GSE144236/norm_log_scale/sptROR2_tcells.rds")
saveRDS(sptnoROR2_macr, file = "./data_output/GSE144236/norm_log_scale/sptnoROR2_macr.rds")
saveRDS(sptnoROR2_tcells, file = "./data_output/GSE144236/norm_log_scale/sptnoROR2_tcells.rds")

## T cells ROR2

sptROR2_tcells@meta.data$sub_clust05 = NA

sptROR2_tcells@meta.data <-
  sptROR2_tcells@meta.data %>%
  mutate(sub_clust05 =
           ifelse(seurat_clusters == "0", "T reg",
                  ifelse(seurat_clusters == "1", "CD44_hi CD8_hi G2MB_hi",
                         ifelse(seurat_clusters == "2", "CD3+",
                                ifelse(seurat_clusters == "3", "CD62L_hi",
                                       ifelse(seurat_clusters == "4", "Monocytes",
                                              NA))))))

Idents(sptROR2_tcells) <- sptROR2_tcells@meta.data$sub_clust05
temp <- DimPlot(sptROR2_tcells, reduction = "tsne")
LabelClusters(temp, id = "ident")

ggsave("tcells_patient_tumor_ror2p_tsne_subclust.png",last_plot(),
       "./figures/GSE144236/", device = png, dpi = 300,
       width = 8, height = 5)

## T cells no ROR2

sptnoROR2_tcells@meta.data$sub_clust05 = NA

sptnoROR2_tcells@meta.data <-
  sptnoROR2_tcells@meta.data %>%
  mutate(sub_clust05 =
           ifelse(seurat_clusters == "0", "CD27_hi TNFR_hi NKG27_hi",
                  ifelse(seurat_clusters == "1", "Natural killer",
                         ifelse(seurat_clusters == "2", "CD4_hi CD44_hi",
                                ifelse(seurat_clusters == "3", "HLA_hi PCNA_hi",
                                              NA)))))

Idents(sptnoROR2_tcells) <- sptnoROR2_tcells@meta.data$sub_clust05
temp <- DimPlot(sptnoROR2_tcells, reduction = "tsne")
LabelClusters(temp, id = "ident")

ggsave("tcells_patient_tumor_ror2n_tsne_subclust.png",last_plot(),
       "./figures/GSE144236/", device = png, dpi = 300,
       width = 8, height = 5)


## Macrophages ROR2

sptROR2_macr@meta.data$sub_clust05 = NA

sptROR2_macr@meta.data <-
  sptROR2_macr@meta.data %>%
  mutate(sub_clust05 =
           ifelse(seurat_clusters == "0", "IL1B_hi LTA4H_hi CLEC12A_hi",
                  ifelse(seurat_clusters == "1", "CD1C_hi HLA_hi",
                         ifelse(seurat_clusters == "2", "B2M_hi CD86_hi",
                                ifelse(seurat_clusters == "3", "LTB_hi CD1A/B_hi CD1A_hi",
                                              NA)))))

Idents(sptROR2_macr) <- sptROR2_macr@meta.data$sub_clust05
temp <- DimPlot(sptROR2_macr, reduction = "tsne")
LabelClusters(temp, id = "ident")

ggsave("macr_patient_tumor_ror2p_tsne_subclust.png",last_plot(),
       "./figures/GSE144236/", device = png, dpi = 300,
       width = 8, height = 5)

## Macrophages no ROR2

sptnoROR2_macr@meta.data$sub_clust05 = NA

sptnoROR2_macr@meta.data <-
  sptnoROR2_macr@meta.data %>%
  mutate(sub_clust05 =
           ifelse(seurat_clusters == "0", "IL1B_hi CD14_hi",
                  ifelse(seurat_clusters == "1", "AIF1_hi B2M_hi",
                         ifelse(seurat_clusters == "2", "CD1C/E_hi HLA_hi S100A8_hi",
                                       NA))))

Idents(sptnoROR2_macr) <- sptnoROR2_macr@meta.data$sub_clust05
temp <- DimPlot(sptnoROR2_macr, reduction = "tsne")
LabelClusters(temp, id = "ident")

ggsave("macr_patient_tumor_ror2n_tsne_subclust.png",last_plot(),
       "./figures/GSE144236/", device = png, dpi = 300,
       width = 8, height = 5)

saveRDS(sptROR2_macr, file = "./data_output/GSE144236/norm_log_scale/sptROR2_macr.rds")
saveRDS(sptROR2_tcells, file = "./data_output/GSE144236/norm_log_scale/sptROR2_tcells.rds")
saveRDS(sptnoROR2_macr, file = "./data_output/GSE144236/norm_log_scale/sptnoROR2_macr.rds")
saveRDS(sptnoROR2_tcells, file = "./data_output/GSE144236/norm_log_scale/sptnoROR2_tcells.rds")

## Getting top10 varying genes in the subclusters of the sROR2_ datasets--------


markers_stcells <- FindAllMarkers(object = sptROR2_tcells,
                          only.pos = F,
                          logfc.threshold = 0.5)
write_csv(markers_stcells,"./data/sptROR2_tcells_markers_diff_per_cluster.csv")

Idents(sptROR2_macr)= "RNA_snn_res.0.5"

top10 =
  markers_stcells %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) 
DoHeatmap(sptROR2_tcells, features = top10$gene) + NoLegend()
write_csv(top10,"./data_output/GSE144236/sptROR2_tcells_markers_diff_per_cluster_top10.csv")
ggsave("./figures/GSE144236/spt_tcells_top10.png", last_plot(), dpi = 300,
       height = 6, width = 9)



markers_stcellsno <- FindAllMarkers(object = sptnoROR2_tcells,
                          only.pos = F,
                          logfc.threshold = 0.5)
write_csv(markers_stcellsno,"./data/sptnoROR2_tcells_markers_diff_per_cluster.csv")
top10 =
  markers_stcellsno %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) 
DoHeatmap(sptnoROR2_tcells, features = top10$gene) + NoLegend()
write_csv(top10,"./data_output/GSE144236/sptnoROR2_tcells_markers_diff_per_cluster_top10.csv")
ggsave("./figures/GSE144236/sptno_tcells_top10.png", last_plot(), dpi = 300,
       height = 6, width = 10)



markers_stmacr <- FindAllMarkers(object = sptROR2_macr,
                          only.pos = F,
                          logfc.threshold = 0.5)
write_csv(markers_stmacr,"./data/sptROR2_macr_markers_diff_per_cluster.csv")
top10 =
  markers_stmacr %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) 
DoHeatmap(sptROR2_macr, features = top10$gene) # + NoLegend()
write_csv(top10,"./data_output/GSE144236/sptROR2_macr_markers_diff_per_cluster_top10.csv")
ggsave("./figures/GSE144236/spt_macr_top10.png", last_plot(), dpi = 300,
       height = 6, width = 9)



markers_stmacrno <- FindAllMarkers(object = sptnoROR2_macr,
                          only.pos = F,
                          logfc.threshold = 0.5)
write_csv(markers_stmacrno,"./data/sptnoROR2_macr_markers_diff_per_cluster.csv")
top10 =
  markers_stmacrno %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) 
DoHeatmap(sptnoROR2_macr, features = top10$gene) + NoLegend()
write_csv(top10,"./data_output/GSE144236/sptnoROR2_macr_markers_diff_per_cluster_top10.csv")
ggsave("./figures/GSE144236/sptno_macr_top10.png", last_plot(), dpi = 300,
       height = 6, width = 9)


## Jingjing meeting
## top10 varying genes, but as dotplots

mnoROR2_macr <- read_csv("data_output/GSE144236/sptnoROR2_macr_markers_diff_per_cluster_top10.csv")
mROR2_macr <- read_csv("data_output/GSE144236/sptROR2_macr_markers_diff_per_cluster_top10.csv")
mnoROR2_tcells <- read_csv("data_output/GSE144236/sptnoROR2_tcells_markers_diff_per_cluster_top10.csv")
mROR2_tcells <- read_csv("data_output/GSE144236/sptROR2_tcells_markers_diff_per_cluster_top10.csv")

# Macrophages top10 dotplot-----------------------------------------------------

mtarget = as.vector(mnoROR2_macr$gene) %>% 
  append(as.vector(mROR2_macr$gene)) %>% 
  append(c("B2M","CD11C","CD11D","ADGRE1","EMR1","IL1B","HLA-DRA","HLA-DRB1","CD274",
           "SOCS1","SOCS3","MRC1","TGM2","CD163","CD40","CD86","CD68",'IDO1'))%>% 
  unique() %>% 
  sort()

dp <- DotPlot(sptnoROR2_macr, features = mtarget, cols = c("red","blue"), dot.scale = 6)+
  FontSize(x.text = 8)
dp + theme(axis.text.x = element_text(angle = 90),
           legend.text=element_text(size=10),
           legend.title=element_text(size=10))+
      labs(title = "Macrophages clusters in epithelial ROR2- patients")
ggsave("./figures/GSE144236/sptno_macr_top10_sorted_dotplot.png", last_plot(), dpi = 300,
       height = 6, width = 11)


# PDL1 (CD274) CD62L (SELL) ADGRE1 (EMR1)

dp <- DotPlot(sptROR2_macr, features = mtarget, cols = c("red","blue"), dot.scale = 6)+
  FontSize(x.text = 8)
dp + theme(axis.text.x = element_text(angle = 90),
           legend.text=element_text(size=10),
           legend.title=element_text(size=10))+
     labs(title = "Macrophages clusters in epithelial ROR2+ patients")
ggsave("./figures/GSE144236/spt_macr_top10_sorted_dotplot.png", last_plot(), dpi = 300,
       height = 6, width = 11)

# Tcells top10 dotplot----------------------------------------------------------

ttarget = as.vector(mnoROR2_tcells$gene) %>% 
  append(as.vector(mROR2_tcells$gene)) %>% 
  append(c("CD8A","CD8B","IL2","HLA-DRA","HLA-DRB1",'IFNG',"GZMB","FOXP3","CD44",
           "CCR7","CCT6A","CD4","CD69",
           "CD62L","SELL","CD274", "XCL1","CD56","CD14","FCGR3A","CD163","CD68","CD86"))%>% 
  unique() %>% 
  sort()

dp <- DotPlot(sptnoROR2_tcells, features = ttarget, cols = c("red","blue"), dot.scale = 6)+
  FontSize(x.text = 8)
dp + theme(axis.text.x = element_text(angle = 90),
           legend.text=element_text(size=10),
           legend.title=element_text(size=10))+
    labs(title = "T cells clusters in epithelial ROR2- patients")
ggsave("./figures/GSE144236/sptno_tcells_top10_dotplot.png", last_plot(), dpi = 300,
       height = 6, width = 11)


dp <- DotPlot(sptROR2_tcells, features = ttarget, cols = c("red","blue"), dot.scale = 6)+
  FontSize(x.text = 8)
dp + theme(axis.text.x = element_text(angle = 90),
           legend.text=element_text(size=10),
           legend.title=element_text(size=10))+
  labs(title = "T cells clusters in epithelial ROR2+ patients")
ggsave("./figures/GSE144236/spt_tcells_top10_dotplot.png", last_plot(), dpi = 300,
       height = 6, width = 11)





























