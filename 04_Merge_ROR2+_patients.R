#####################################
## Merge all ROR2+ tumor patients ##
####################################

library(Seurat)
library(tidyverse)

dir <- getwd()
ldir <- list.files(paste0(dir,"/data_output/GSE144236/norm_log_scale/"), pattern = "seurat_f.rds")
lnames <- c("P1_Normal","P1_Tumor","P10_Normal","P10_Tumor","P2_Normal","P2_Tumor",
            "P3_Normal","P3_Tumor","P4_Normal","P4_Tumor","P5_Normal","P5_Tumor",
            "P6_Normal","P6_Tumor","P7_Normal","P7_Tumor","P8_Normal",
            "P8_Tumor","P9_Normal","P9_Tumor")

for(i in 1:length(ldir)){
  temp <- readRDS(paste0(dir,"/data_output/GSE144236/norm_log_scale/", ldir[i]))
  temp@meta.data$patient <- lnames[i]
  temp@meta.data$condition <- sub(".*_", "",lnames[i])
  assign(paste0(lnames[i],"_seurat_f"), temp)
}

seurat_all <- merge(P1_Tumor_seurat_f, c(P1_Normal_seurat_f, P10_Normal_seurat_f,
                                        P10_Tumor_seurat_f,P2_Normal_seurat_f,P2_Tumor_seurat_f,
                                        P3_Tumor_seurat_f,P3_Normal_seurat_f,P4_Normal_seurat_f,
                                        P4_Tumor_seurat_f,P5_Normal_seurat_f,P5_Tumor_seurat_f,
                                        P6_Normal_seurat_f,P6_Tumor_seurat_f,P7_Normal_seurat_f,
                                        P7_Tumor_seurat_f,P8_Normal_seurat_f,
                                        P8_Tumor_seurat_f,P9_Normal_seurat_f,P9_Tumor_seurat_f))


table(seurat_all@meta.data$patient)
length(P1_Tumor_seurat_f@meta.data$patient)

saveRDS(seurat_all,"./data_output/GSE144236/norm_log_scale/seurat_all.rds")

rm(list = c("P1_Normal_seurat_f","P10_Normal_seurat_f",
            "P10_Tumor_seurat_f","P2_Normal_seurat_f","P2_Tumor_seurat_f",
            "P3_Normal_seurat_f","P3_Tumor_seurat_f","P4_Normal_seurat_f",
            "P4_Tumor_seurat_f","P5_Normal_seurat_f","P5_Tumor_seurat_f",
            "P6_Normal_seurat_f","P6_Tumor_seurat_f","P7_Normal_seurat_f",
            "P7_Tumor_seurat_f","P8_Normal_seurat_f",
            "P8_Tumor_seurat_f","P9_Normal_seurat_f","P9_Tumor_seurat_f", "temp",
            "P1_Normal_seurat_f","P1_Tumor_seurat_f"))

# Clustering

seurat_all <- FindVariableFeatures(seurat_all, selection.method = "vst", nfeatures = 2000)

s.genes <- read_csv("./data/GSE144236/S_phase_gene.csv")
g2m.genes <- read.csv("./data/GSE144236/G2M_phase_gene.csv")

seurat_all <- ScaleData(seurat_all, vars.to.regress = 
                        c("S.Score", "G2M.Score","mitoPercent","nCount_RNA"), 
                      features = rownames(seurat_all))

seurat_all <- RunPCA(seurat_all, features = VariableFeatures(object = seurat_all))

seurat_all <- FindNeighbors(seurat_all, dims = 1:14, verbose = FALSE)
seurat_all <- FindClusters(seurat_all, resolution = 0.5, verbose = FALSE)

seurat_all <- RunUMAP(seurat_all, dims = 1:20, verbose = FALSE)
DimPlot(seurat_all, reduction = "umap", label = T)
Idents(seurat_all) = "seurat_clusters"

plot <- DimPlot(object = seurat_all)
LabelClusters(plot = plot, id = 'ident')
ggsave("ALL_UMAP_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png", height = 7, width = 10)

seurat_all <- RunTSNE(seurat_all, dims = 1:20, verbose = FALSE)
DimPlot(seurat_all, reduction = "tsne")
plot1 <- DimPlot(seurat_all, reduction = "tsne")
LabelClusters(plot1, id = "ident")

ggsave("ALL_tsne_labelled.png",last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png", height = 7, width = 10)

# https://distill.pub/2016/misread-tsne/

saveRDS(seurat_all, file = "./data_output/GSE144236/norm_log_scale/seurat_all_processed.rds")

p1 <- DimPlot(seurat_all, reduction = "tsne", group.by = "patient")
p2 <- DimPlot(seurat_all, reduction = "tsne", group.by = "condition", 
              label = TRUE,
              repel = TRUE)
p1 + p2

ggsave("ALL_tsne_patient_condition.png",last_plot(),"./figures/GSE144236/", dpi = 500 ,
       device = "png", width = 22, height = 10)




