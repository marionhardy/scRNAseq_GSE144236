#######################################################
## Reprocess individual data but using scale and log ##
#######################################################

library(Seurat)
library(tidyverse)


lnames <- c("P1_Normal","P1_Tumor","P10_Normal","P10_Tumor","P2_Normal","P2_Tumor",
            "P3_Normal","P3_Tumor","P4_Normal","P4_Tumor","P5_Normal","P5_Tumor",
            "P6_Normal","P6_Tumor","P7_Normal","P7_Tumor","P8_Normal",
            "P8_Tumor","P9_Normal","P9_Tumor")


s.genes <- read_csv("./data/GSE144236/S_phase_gene.csv")
g2m.genes <- read.csv("./data/GSE144236/G2M_phase_gene.csv")

# Make the Seurat object with the raw data

dir <- getwd()
ldir <- list.files(paste0(dir,"/data_output/GSE144236/"), pattern = ".csv")

for(i in 1:length(ldir)){
  
  p <- read.csv(paste0(dir,"/data_output/GSE144236/", ldir[i]), row.names = 1)
  
  seurat <- CreateSeuratObject(counts = p, project = "221122", min.cells = 3, 
                               min.features = 10)
  seurat@meta.data$patient <- lnames[i]
  seurat@meta.data$condition <- sub(".*_", "",lnames[i])
  setwd(dir)
  
  seurat$mitoPercent <- PercentageFeatureSet(seurat, pattern='^MT-')

  seurat_f <- subset(seurat, subset = nFeature_RNA > 200 & 
                       nFeature_RNA < 6000 & 
                       mitoPercent < 10)
  rm(list = "seurat")
  rm(list = "p")
  seurat_f <- NormalizeData(seurat_f, normalization.method = "LogNormalize")
  seurat_f <- FindVariableFeatures(seurat_f, selection.method = "vst", nfeatures = 2000)
  seurat_f <- ScaleData(seurat_f, features = rownames(seurat_f))

  assign(paste0(lnames[i],"_seurat_f"),seurat_f)
  saveRDS(seurat_f, file = paste0("./data_output/GSE144236/norm_log_scale/",
                                  lnames[i],"_seurat_f.rds"))
}






















