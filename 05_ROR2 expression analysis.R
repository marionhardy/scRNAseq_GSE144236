
library(Seurat)
library(tidyverse)
library(RColorBrewer)
library(gridExtra)
library(viridis)
library(Polychrome)
library(circlize)
library(ComplexHeatmap)
library(patchwork)

seurat_all <-  read_rds(file = "./data_output/GSE144236/norm_log_scale/seurat_all_processed.rds")

# Find markers for each clusters------------------------------------------------

Idents(seurat_all) = "RNA_snn_res.0.5"
markers <- FindAllMarkers(object = seurat_all,
                          only.pos = F,
                          logfc.threshold = 0.5)

write_csv(markers,"./data_output/GSE144236/norm_log_scale/ALL_markers_diff_per_cluster.csv")

## Juliette and Chim-Kei's questions--------------------------------------------

target <- c("ROR2","ROR1","STAT1","STAT3","CTNNB1")

FeaturePlot(object = seurat_all,
            features = "ROR2",
            reduction = "tsne",
            sort.cell = TRUE,
            min.cutoff = 'q10', 
            label = TRUE,
            repel = TRUE) # + 
  #scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "YlOrRd")))


ggsave("ALL_ROR2_tsne.png", plot = last_plot(), path = "./figures/GSE144236/",
       device = "png", dpi = 500, width = 5, height = 4)

VlnPlot(object = seurat_all, features = "ROR2")

ggsave("ALL_ROR2_violinplot.png", plot = last_plot(), path = "./figures/GSE144236/",
       device = "png", dpi = 500, width = 10, height = 5)


# Markers for celltype identification

auth_celltype <- c("KRT1","KRT5","KRT10","KRT14","LYZ","HLA-DRB1", "HLA-DRA","HLA-DQB2",
              "CD3D", "CD2","CD7","COL1A1", "COL1A2","LUM","MLANA", "DCT", "PMEL",
              "TFF3", "CLDN5", "VWF", "IGLL5", "IGJ"," MS4A1", "CD79A","CD163",
              "CD68","S100A8"," S100A9", "TREM1","CD1C", "CLEC10A","CLEC9A", "CADM1", 
              "XCR1", "CD207", "CD1A", "S100B","AXL", "IGFBP5", "PPP1R14A","CLEC4C",
              "IL3RA","CCR7", "CCL19","CD4","CD8")

celltype <- c("KRT1","KRT5","KRT10","KRT14","LYZ","HLA-DRB1", "HLA-DRA","HLA-DQB2",
              "HLA-DQB1","B2M","CD1C","CD11C","CD11B","CLEC9A","XCR1", "CD207",
              "S100B","S100A9","ADGRE1","EMR1","CD68","CD86","CD163","TFF3","CLDN5","VWF",
              "MS4A1","CD79A","IGJ","IGLL5","CD3D", "CD2","CD7", "CD3G","CD4","CD8A",
              "CD8B","CD56","NCAM","CXCR1","CXCR3","XCL1","COL1A1","COL1A2","LUM",
              "MLANA","PMEL","DCT","CD14","FCGR3A")

Idents(seurat_all) = "RNA_snn_res.0.5"
dp <- DotPlot(seurat_all, features = celltype, cols = c("red","blue"), dot.scale = 6)+
  FontSize(x.text = 8)
dp + theme(axis.text.x = element_text(angle = 90),
           legend.text=element_text(size=10),
           legend.title=element_text(size=10))+
  labs(title = "Cell specific clustering markers for seurat_all")
ggsave("./figures/GSE144236/ALL_clustering_markers.png", last_plot(), dpi = 300,
       height = 6, width = 11)

#  Visualize with complexheatmap

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

png(filename ="./figures/GSE144236/ALL_heatmap.png")

plot <- Heatmap(exp_mat,
                heatmap_legend_param=list(title="expression"),
                column_title = "clustered dotplot", 
                col=col_fun,
                rect_gp = gpar(type = "none"),
                cell_fun = cell_fun,
                row_names_gp = gpar(fontsize = 6),
                row_km = 5,
                border = "black")
print(plot)
dev.off() 

## Attribute cell types to clusters

seurat_all@meta.data$attr_clusters <- NA
seurat_all@meta.data$attr_gr_clusters <- NA
seurat_all@meta.data$paper_clusters <- NA


seurat_all@meta.data <-
  seurat_all@meta.data %>%
  mutate(attr_clusters=
           ifelse(seurat_clusters == "0", "Epithelial_1",
                  ifelse(seurat_clusters == "1", "Epithelial_2",
                         ifelse(seurat_clusters == "2", "Epithelial_3",
                                ifelse(seurat_clusters == "3", "DCs_1",
                                       ifelse(seurat_clusters == "4", "Epithelial_4",
                                              ifelse(seurat_clusters == "5", "Macrophages_1",
                                                     ifelse(seurat_clusters == "6", "Langerhan cells_1", 
                                                            ifelse(seurat_clusters == "7", "Epithelial_5",
                                                                   ifelse(seurat_clusters == "8", "Macrophages_2",
                                                                          ifelse(seurat_clusters == "9", "T cells",
                                                                                 ifelse(seurat_clusters == "10", "MDSCs",
                                                                                        ifelse(seurat_clusters == "11", "Langerhan cells_2",
                                                                                               ifelse(seurat_clusters == "12", "Melanocytes",
                                                                                                      ifelse(seurat_clusters == "13", "Fibroblasts_1",
                                                                                                             ifelse(seurat_clusters == "14", "DCs_2", 
                                                                                                                    ifelse(seurat_clusters == "15", "Endothelial",
                                                                                                                           ifelse(seurat_clusters == "16", "Epithelial_6",
                                                                                                                                  ifelse(seurat_clusters == "17", "B cell/plasma",
                                                                                                                                         ifelse(seurat_clusters == "18", "B cell/plasma",
                                                                                                                                                ifelse(seurat_clusters == "19", "Fibroblasts_2",
                                                                                                                                                NA)))))))))))))))))))))

Idents(seurat_all) = "attr_clusters"
DimPlot(seurat_all, reduction = "tsne", label = TRUE, pt.size = 0.3)

ggsave("ALL_tsne_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png", width = 14, height = 10)


seurat_all@meta.data <-
  seurat_all@meta.data %>%
  mutate(attr_gr_clusters=
           ifelse(seurat_clusters %in% c("0","1","2","4","7","16"), "Epithelial cells",
                  ifelse(seurat_clusters %in% c("5","8"), "Macrophages",
                         ifelse(seurat_clusters == "9", "T cells",
                                ifelse(seurat_clusters %in% c("13","19"), "Fibroblasts",
                                       ifelse(seurat_clusters == "12", "Melanocytes",
                                                     ifelse(seurat_clusters == "15", "Endothelial cells", 
                                                            ifelse(seurat_clusters %in% c("17","18"), "B cell/plasma",
                                                                   ifelse(seurat_clusters %in% c("14","3"), "DCs",
                                                                          ifelse(seurat_clusters %in% c("10"), "MDSCs",
                                                                                 ifelse(seurat_clusters %in% c("6","11"), "Langerhan cells",
                                                                   NA)))))))))))

Idents(seurat_all) <- seurat_all@meta.data$attr_gr_clusters
temp <- DimPlot(seurat_all, reduction = "tsne")
LabelClusters(temp, id = "ident")

ggsave("ALL_tsne_labelled_grouped_attr.png",last_plot(),"./figures/GSE144236/", device = png, dpi = 500,
       width = 11, height = 8)


seurat_all@meta.data <-
  seurat_all@meta.data %>%
  mutate(paper_clusters=
           ifelse(seurat_clusters %in% c("0","1","2","4","7","16"), "Epithelial cells",
                  ifelse(seurat_clusters %in% c("3","14","5","8","10","6","11"), "Myeloid cells",
                         ifelse(seurat_clusters == "9", "T cells",
                                ifelse(seurat_clusters %in% c("13","19"), "Fibroblasts",
                                       ifelse(seurat_clusters == "12", "Melanocytes",
                                              ifelse(seurat_clusters == "15", "Endothelial cells", 
                                                     ifelse(seurat_clusters %in% c("17","18"), "B cell/plasma",
                                                            NA))))))))

Idents(seurat_all) <- seurat_all@meta.data$paper_clusters
temp <- DimPlot(seurat_all, reduction = "tsne")
LabelClusters(temp, id = "ident")

ggsave("ALL_tsne_labelled_paper.png",last_plot(),"./figures/GSE144236/", device = png, dpi = 500,
       width = 11, height = 8)

saveRDS(seurat_all, file = "./data_output/GSE144236/norm_log_scale/seurat_all_processed.rds")

# Get some stats babyyy

## Global

table(seurat_all@meta.data$paper_clusters, seurat_all@meta.data$patient)
table(seurat_all@meta.data$attr_clusters, seurat_all@meta.data$patient)
table(seurat_all@meta.data$patient, seurat_all@meta.data$condition)

# nbr of cells per cluster
table(Idents(seurat_all))
# nbr of ROR2+ cells
table(Idents(subset(seurat_all, ROR2 > 0))) 

# average expression level of ROR2 per patient
VlnPlot(object = seurat_all, features = "ROR2")
VlnPlot(object = seurat_all, features = "ROR2", split.by = "condition")

# Target now with labelled clusters

table(seurat_all@meta.data$paper_clusters)

FeaturePlot(object = seurat_all,
            reduction = "tsne",
            features = "ROR2",
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 4,
            repel = TRUE,
            keep.scale = "all")

ggsave("ALL_ROR2_tSNE_labelled.png",last_plot(),"./figures/GSE144236/", device = png, dpi = 500,
       width = 5, height = 4)

VlnPlot(object = seurat_all, features = "ROR2")
ggsave("ALL_ROR2_violinplots_labelled.png",last_plot(),"./figures/GSE144236/", device = png, 
       dpi = 500)

VlnPlot(object = seurat_all, features = "ROR2", split.by = "condition", fill.by = "ident")
ggsave("ALL_ROR2_violinplots_labelled.png",last_plot(),"./figures/GSE144236/", device = png, 
       dpi = 500)


## Subset on cells that are ROR2+

sROR2 <- subset(seurat_all, subset = ROR2 > 0)

saveRDS(sROR2, file = "./data_output/GSE144236/norm_log_scale/seurat_all_ROR2pos.rds")

VlnPlot(object = sROR2, features = "ROR2", split.by = "condition", fill.by = "ident")
ggsave("ALL_ROR2pos_violinplots_labelled.png",last_plot(),"./figures/GSE144236/", device = png, 
       dpi = 500, width = 8, height = 6)

VlnPlot(object = sROR2, features = "ROR2", split.by = "patient", fill.by = "ident")
ggsave("ALL_ROR2pos_violinplots_patients.png",last_plot(),"./figures/GSE144236/", device = png, 
       dpi = 500, width = 12, height = 6)

## Subset per Normal or tumor 

DimPlot(object = seurat_all, split.by = "condition", reduction = "umap")

ggsave("SPLIT_umap_labelled_paper.png",last_plot(),"./figures/GSE144236/", device = png, dpi = 500,
       width = 16, height = 8)


DimPlot(object = seurat_all, split.by = "condition", reduction = "tsne")

ggsave("SPLIT_tsne_labelled_paper.png",last_plot(),"./figures/GSE144236/", device = png, dpi = 500,
       width = 16, height = 8)

FeaturePlot(object = seurat_all,
            reduction = "tsne",
            features = "ROR2",
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3,
            pt.size = 0.6,
            repel = TRUE,
            keep.scale = "all",
            split.by = "condition")

ggsave("SPLIT_ROR2_tSNE_labelled.png",last_plot(),"./figures/GSE144236/", device = png, dpi = 500,
       width = 14, height = 6)

FeaturePlot(object = seurat_all,
            reduction = "umap",
            features = "ROR2",
            sort.cell = TRUE,
            min.cutoff = 'q10',
            label = TRUE,
            label.size = 3.5,
            pt.size = 0.8,
            repel = TRUE,
            keep.scale = "all",
            split.by = "condition")

ggsave("SPLIT_ROR2_UMAP_labelled.png",last_plot(),"./figures/GSE144236/", device = png, dpi = 500,
       width = 14, height = 6)

# Add annotation TRUE/FALSE for cell expressing ROR2 or not
# violin plot of expression level
# ranking of tumors expression levels

sROR2 <- subset(seurat_all, subset = ROR2 > 0)
seurat_all@meta.data$exprROR2 <- ifelse(rownames(seurat_all@meta.data) %in% colnames(sROR2),"Expr", "None")
table(seurat_all@meta.data$exprROR2)
length(colnames(sROR2))


DimPlot(object = seurat_all, split.by = "exprROR2", reduction = "tsne")

ggsave("SPLIT_ROR2_tsne_labelled_paper.png",last_plot(),"./figures/GSE144236/", 
       device = png, dpi = 500,
       width = 16, height = 8)



sROR2_e <- subset(sROR2, paper_clusters=="Epithelial cells")

p <- VlnPlot(object = sROR2_e, features = "ROR2", split.by = "condition",
        group.by= "patient",fill.by = "ident")+
  geom_boxplot(width=0.3, fill="white")

print(p)

ggsave("Epithelial_patient_ROR2_violinboxplots.png",p,"./figures/GSE144236/", device = png, 
       dpi = 500, width = 12, height = 6)


boxp <- as.data.frame(p$data)
boxp %>% 
  ggplot(aes(x=ident, y=ROR2, fill=split)) +
  geom_boxplot()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave("Epithelial_patient_ROR2_boxplots.png",last_plot(),"./figures/GSE144236/", device = png, 
       dpi = 500, width = 12, height = 6)
boxp %>% 
  ggplot(aes(x=ident, y=ROR2, fill=split)) +
  geom_boxplot()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  stat_summary(fun = median, geom = "text", col = "black",     # Add text to plot
               vjust = 1.5, aes(label = paste(round(..y.., digits = 3))))

ggsave("Epithelial_patient_ROR2_labelled_boxplots.png",last_plot(),"./figures/GSE144236/", device = png, 
       dpi = 500, width = 12, height = 6)


# Rank patients based on ROR2 expr level in epithelial tumor cells
# Check activation markers

seurat_all@meta.data$rtpatient <- NA
seurat_all@meta.data$rtpatient <- seurat_all@meta.data$orig.ident
seurat_all@meta.data$rtpatient <- seurat_all@meta.data$patient

seurat_all$rtpatient <- 
  factor(seurat_all$rtpatient, levels = c("P3_Normal","P3_Tumor","P8_Tumor",
                                          "P5_Tumor","P9_Tumor", "P1_Tumor",
                                          "P6_Tumor","P7_Tumor","P2_Tumor",
                                          "P10_Tumor","P5_Normal","P4_Tumor","P10_Normal",
                                          "P9_Normal","P4_Normal","P1_Normal","P7_Normal","P2_Normal",
                                          "P6_Normal","P8_Normal"))

sROR2 <- subset(seurat_all, subset = ROR2 > 0)
snoROR2 <- subset(seurat_all, subset = ROR2 == 0)
sROR2_e <- subset(sROR2, paper_clusters=="Epithelial cells")
sROR2_m <- subset(sROR2, paper_clusters=="Myeloid cells")
seurat_mpos <- subset(seurat_all, paper_clusters=="Myeloid cells"&exprROR2=="Expr")
seurat_mneg <- subset(seurat_all, paper_clusters=="Myeloid cells"&exprROR2=="None")

act_marker <- c("ROR2","CCR7","FOXP3","CD4","CD8A","MKI67","XCL1","CD44","LY6C",
                "CCT6A","CD45RA","CD45RO")
Idents(seurat_all) <- seurat_all$rtpatient

dot1 <- 
  DotPlot(sROR2_e, features = "ROR2", group.by = "rtpatient") +
  scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "firebrick3")+
  NoLegend()+
  labs(title = "Epithelial cells", subtitle = "ROR2+")# Epithelial ROR2 ordered

Idents(seurat_mneg) <- seurat_mneg@meta.data$rtpatient
Idents(seurat_mpos) <- seurat_mpos@meta.data$rtpatient

dot2 <- 
DotPlot(seurat_mpos, features = act_marker, group.by = "rtpatient")+ 
  scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "firebrick3") +
  labs(title = "Myeloid cells", subtitle = "ROR2+")# Activation marker for myeloid cells

dot3 <- 
  DotPlot(seurat_mneg, features = act_marker, group.by = "rtpatient")+ 
  scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "firebrick3") +
  NoLegend()+
  labs(title = "Myeloid cells", subtitle = "ROR2-")# Activation marker for myeloid cells

dot1/dot2/dot3

ggsave("SPLIT_Activation_marker_Dotplots.png", last_plot(),device = "png", dpi = 300,
       height = 4, width = 10, path = "./figures/GSE144236/")



