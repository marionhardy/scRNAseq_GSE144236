#################################
## Chim-Kei and Juliette scRNA ##
#################################

# Based on GSE144240 : Multimodal analysis and spatial architecture in human 
# squamous cell carcinoma

# GSE144236 contains the scRNAseq data only
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA603103&o=acc_s%3Aa
# Raw files

library("GEOquery")
library(tidyverse)

gse <-  getGEO("GSE144236", GSEMatrix = F)
head(Meta(gse))
names(GSMList(gse))
GSMList(gse)[[1]]
names(GPLList(gse))
GPLList(gse)[[1]]

gse144236 <- getGEO('GSE144236',GSEMatrix=TRUE)
show(gse144236)
class(gse144236)
show(pData(phenoData(gse144236[[1]]))[1:28,c(1,6,8)])

# GSM contains all the samples I want to analyze
# However, I don't need the normal samples rn

gsm <- gse144236[[1]]
exprs(gsm) # We are missing the count data
pData(gsm) # It's available as a .tsv supplementary file for each GSMxxxxxxx


# How to extract everything in the same place?

supp_files <- getGEOSuppFiles(GEO = "GSE144236", 
                              baseDir = paste0(getwd(),"/data/")) 
# or download them manually

gunzip("./data/GSE144236/GSE144236_cSCC_counts.txt.gz", overwrite = F, remove = FALSE)

counts <- read.table("./data/GSE144236/GSE144236_cSCC_counts.txt", row.names = 1)

head(counts)[1:5,1:5]

# Extract scRNA counts for each patient

lnames <- c("P1_Tumor","P1_Normal","P2_Tumor","P2_Normal","P3_Tumor",
            "P3_Normal","P4_Tumor","P4_Normal","P5_Tumor","P5_Normal",
            "P6_Tumor","P6_Normal","P7_Normal","P7_Tumor","P8_Tumor",
            "P8_Normal","P9_Tumor","P9_Normal","P10_Tumor","P10_Normal")


for (i in 1:length(lnames)){
  p <- counts %>% select(matches(paste0(lnames[i])))
  p <- p[-c(1:2),]
  write.csv(p,paste0("./data_output/GSE144236/",lnames[i],".csv"))
}

## GSE file for SCC13 xenograft-------------------------------------------------

gunzip("./data/GSE144236/GSE144236_SCC13_counts.txt.gz", overwrite = F, remove = FALSE)
counts <- read.table("./data/GSE144236/GSE144236_SCC13_counts.txt", row.names = 1)

head(counts)[1:5,1:5]

write.csv(counts,"./data_output/GSE144236_SCC13/SCC13_counts.csv")


## GSE file for CAL27 xenograft-------------------------------------------------

gunzip("./data/GSE144236/GSE144236_CAL27_counts.txt.gz", overwrite = F, remove = FALSE)
counts <- read.table("./data/GSE144236/GSE144236_CAL27_counts.txt", row.names = 1)

head(counts)[1:5,1:5]

write.csv(counts,"./data_output/GSE144236_CAL27/CAL27_counts.csv")















































