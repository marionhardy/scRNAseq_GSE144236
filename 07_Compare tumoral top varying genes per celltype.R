
library(Seurat)
library(tidyverse)

# Create seurat with tumoral data only

spt = subset(seurat_all, subset = condition == "Tumor")

# First add metadata ROR2+ and ROR2-

spt@meta.data$exprROR2 <- ifelse(rownames(spt@meta.data) %in% 
                                          colnames(sptROR2),"ROR2+", "ROR2-")
table(spt@meta.data$exprROR2)
# saveRDS(spt,"./data_output/GSE144236/norm_log_scale/spt.rds")

# Compare P3+P5+P8 (they are ROR2-) vs all other patients per cluster

Idents(spt) = spt@meta.data$seurat_clusters

spt$celltype.stim <- paste(Idents(spt), spt$exprROR2, sep = "_")
spt$celltype <- Idents(spt)
Idents(spt) <- "celltype.stim"

for (i in 0:18){ #or however many clusters you have
  try({
    ident1 <- paste0(i,"_ROR2+")
    ident2 <- paste0(i,"_ROR2-")
    condition_diffgenes <- FindMarkers(spt, ident.1 = ident1, ident.2=ident2, 
                                       min.pct=0.25, logfc.threshold=0)
  })
}

condition_diffgenes$genes = rownames(condition_diffgenes)

write_csv(condition_diffgenes,"./data_output/GSE144236/spt_markers_diff_per_attr_cluster.csv")


condition.top10 =
  condition_diffgenes %>% filter(p_val_adj<0.05 & avg_log2FC<=-0.5 |avg_log2FC>=0.5)
target = as.vector(rownames(condition.top10))

FeaturePlot(spt, features = target, split.by = "exprROR2", max.cutoff = 3, 
            cols = c("grey", "red"))
ggsave("./figures/GSE144236/spt_top13_UMAP.png", last_plot(), dpi = 600,
       height = 45, width = 12)


plots = VlnPlot(spt, features = target, split.by = "exprROR2", group.by = "attr_clusters", 
                 pt.size = 0, combine = FALSE, split.plot = T)

wrap_plots(plots = plots, ncol = 1)

ggsave("./figures/GSE144236/spt_top13_violinplot.png", last_plot(), dpi = 600,
       height = 45, width = 12)


# GSEA analysis of varying genes between ROR2+ and ROR2- patient tumors

# library("org.Mm.eg.db") ! we're in a human organism
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")

# GSEA--------------------------------------------------------------------------

condition_diffgenes_f = condition_diffgenes %>% 
  filter(p_val_adj <= 0.05)
ordered_genes <- abs(condition_diffgenes_f$avg_log2FC)
names(ordered_genes) <- condition_diffgenes_f$genes
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

# GO

go_gsea <- gseGO(gene = ordered_genes,
                 OrgDb = org.Hs.eg.db,
                 scoreType = "pos",
                 keyType = "SYMBOL",
                 ont          = "ALL",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

dotplot(go_gsea, showCategory=30) + ggtitle("Dotplot for GSEA all")
ggsave( "./figures/GSE144236/GSEA_ALL.png", plot = last_plot(), device = "png", dpi = 300,
        width = 8, height = 7)

dotplot(go_gsea %>% filter(ONTOLOGY=="BP"), showCategory=30)+ ggtitle("Dotplot for GSEA BP")
ggsave("./figures/GSE144236/GSEA_BP.png", plot = last_plot(), device = png, dpi = 300,
       width = 6, height = 9)

dotplot(go_gsea, split = "ONTOLOGY",showCategory=10, x="NES") + 
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GSEA all split")
ggsave("./figures/GSE144236/GSEA_ALL_split.png", plot = last_plot(), device = png, dpi = 300,
       width = 8, height = 12)

go_gsea_tbl <- as.tibble(go_gsea)


# ------------------------------------------------------------------------------

## If we filter and only compare the macrophages, the epithelial cells and the 
# t cells ROR2+ vs ROR2-?

Idents(spt) = spt@meta.data$paper_clusters

spt$celltype.stim.auth <- paste(Idents(spt), spt$exprROR2, sep = "_")
spt$celltype <- Idents(spt)
Idents(spt) <- "celltype.stim.auth"
table(spt@meta.data$celltype.stim.auth)

# By hand

cond.diff.epith <- FindMarkers(spt, ident.1 = "Epithelial cells_ROR2+", 
                              ident.2="Epithelial cells_ROR2-",
                              min.pct=0.25, logfc.threshold=0)

cond.diff.epith$genes = rownames(cond.diff.epith)

write_csv(cond.diff.epith,"./data_output/GSE144236/norm_log_scale/spt_markers_diff_per_epithelial_cells.csv")


# GSEA 

# library("org.Mm.eg.db") ! we're in a human organism
library("org.Hs.eg.db")
library("clusterProfiler")
library("enrichplot")

# GSEA--------------------------------------------------------------------------

cond.diff.epith_f = cond.diff.epith %>% 
  filter(p_val_adj <= 0.05)
ordered_genes <- abs(cond.diff.epith_f$avg_log2FC)
names(ordered_genes) <- cond.diff.epith_f$genes
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

# GO

go_gsea <- gseGO(gene = ordered_genes,
                 OrgDb = org.Hs.eg.db,
                 scoreType = "pos",
                 keyType = "SYMBOL",
                 ont          = "ALL",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

dotplot(go_gsea, showCategory=30) + ggtitle("Dotplot for GSEA all epithelial")
ggsave( "./figures/GSE144236/GSEA_epith_ALL.png", plot = last_plot(), device = "png", dpi = 300,
        width = 8, height = 7)

dotplot(go_gsea %>% filter(ONTOLOGY=="BP"), showCategory=30)+ ggtitle("Dotplot for GSEA BP epithelial")
ggsave("./figures/GSE144236/GSEA_epith_BP.png", plot = last_plot(), device = png, dpi = 300,
       width = 6, height = 9)

go_gsea_tbl <- as.tibble(go_gsea)


## For myeloid cells

# By hand

cond.diff.macr <- FindMarkers(spt, ident.1 = "Myeloid cells_ROR2+", 
                               ident.2="Myeloid cells_ROR2-",
                               min.pct=0.25, logfc.threshold=0)

cond.diff.macr$genes = rownames(cond.diff.macr)

write_csv(cond.diff.macr,"./data_output/GSE144236/norm_log_scale/spt_markers_diff_per_myeloid.csv")


cond.diff.macr_f = cond.diff.macr %>% 
  filter(p_val_adj <= 0.05)
ordered_genes <- abs(cond.diff.macr_f$avg_log2FC)
names(ordered_genes) <- cond.diff.macr_f$genes
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

# GO

go_gsea <- gseGO(gene = ordered_genes,
                 OrgDb = org.Hs.eg.db,
                 scoreType = "pos",
                 keyType = "SYMBOL",
                 ont          = "ALL",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)

dotplot(go_gsea, showCategory=30) + ggtitle("Dotplot for GSEA all myeloid cells")
ggsave( "./figures/GSE144236/myeloid/GSEA_myeloid_ALL.png", plot = last_plot(), device = "png", dpi = 300,
        width = 8, height = 7)

dotplot(go_gsea %>% filter(ONTOLOGY=="BP"), showCategory=30, x="NES")+ ggtitle("Dotplot for GSEA BP myeloid cells")
ggsave("./figures/GSE144236/myeloid/GSEA_myeloid_BP.png", plot = last_plot(), device = png, dpi = 300,
       width = 6, height = 9)

dotplot(go_gsea, split = "ONTOLOGY",showCategory=10, x="NES") + 
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GSEA all split myeloid cells")
ggsave("./figures/GSE144236/myeloid/GSEA_myeloid_ALL_split.png", plot = last_plot(), device = png, dpi = 300,
       width = 8, height = 12)

go_gsea_tbl <- as.tibble(go_gsea)

# MsigDb

library(msigdbr)

hsa_reactome_sets <- msigdbr(
  species = "Homo sapiens", 
  category = "C2",
  subcategory = "CP:REACTOME") # for reactome collection

hsa_kegg_sets <- msigdbr(
  species = "Homo sapiens", 
  category = "C2",
  subcategory = "CP:KEGG") # for KEGG collection

hsa_wiki_sets <- msigdbr(
  species = "Homo sapiens", 
  category = "C2",
  subcategory = "CP:WIKIPATHWAYS") # for Wikipathways collection

set.seed(69)

gsea_results_react <- GSEA(
  geneList = ordered_genes, 
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(
    hsa_reactome_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
) 

dotplot(gsea_results_react, x = "NES", showCategory = 30)+ ggtitle("GSEA reactome LFC myeloid cells")
ggsave("./figures/GSE144236/myeloid/dotplot_GSEA_reactome_myeloid.png", plot = last_plot(), device = png, dpi = 400,
       width = 10, height = 8)

gsea_results_kegg <- GSEA(
  geneList = ordered_genes,
  pvalueCutoff = 0.05, 
  eps = 0,
  seed = TRUE, 
  pAdjustMethod = "BH",
  scoreType = "pos",
  TERM2GENE = dplyr::select(
    hsa_kegg_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_kegg, x = "NES", showCategory = 30)+ ggtitle("GSEA reactome LFC myeloid cells")
ggsave("./figures/GSE144236/myeloid/dotplot_GSEA_kegg_myeloid.png", plot = last_plot(), device = png, dpi = 400,
       width = 10, height = 8)

# It's still possible to use pathview to extract the significant paths

library("pathview")

keggresids <- c("04060","04612")

foldchanges <- cond.diff.macr_f$avg_log2FC
names(foldchanges) <- cond.diff.macr_f$genes
head(foldchanges)
table(is.na(foldchanges))

tmp <- sapply(keggresids, function(pid) pathview(gene.data = foldchanges,
                                                 gene.idtype = "SYMBOL",
                                                 pathway.id = pid,
                                                 species = "hsa",
                                                 kegg.dir="./figures/GSE144236/myeloid",
                                                 out.suffix= "_Colored",
                                                 kegg.native = TRUE,
                                                 map.null = FALSE))


gsea_results_wiki <- GSEA(
  geneList = ordered_genes, 
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    hsa_wiki_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_wiki, x = "NES", showCategory = 30)+ ggtitle("GSEA Wikipathway LFC myeloid cells")

gsea_results_wiki <- as.tibble(gsea_results_wiki)

ggsave("./figures/GSE144236/myeloid/dotplot_GSEA_Wikipathway.png", plot = last_plot(), device = png, dpi = 400)


## For T cells

# By hand

cond.diff.tcell <- FindMarkers(spt, ident.1 = "T cells_ROR2+", 
                              ident.2="T cells_ROR2-",
                              min.pct=0.25, logfc.threshold=0)

cond.diff.tcell$genes = rownames(cond.diff.tcell)

write_csv(cond.diff.tcell,"./data_output/GSE144236/norm_log_scale/spt_markers_diff_per_tcells.csv")


cond.diff.tcell_f = cond.diff.tcell %>% 
  filter(p_val_adj <= 0.05)
ordered_genes <- abs(cond.diff.tcell_f$avg_log2FC)
names(ordered_genes) <- cond.diff.tcell_f$genes
ordered_genes <- sort(ordered_genes, decreasing = TRUE)

# GO

go_gsea <- gseGO(gene = ordered_genes,
                 OrgDb = org.Hs.eg.db,
                 scoreType = "pos",
                 keyType = "SYMBOL",
                 ont          = "ALL",
                 minGSSize    = 100,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = T)

dotplot(go_gsea, showCategory=30) + ggtitle("Dotplot for GSEA all")
ggsave( "./figures/GSE144236/tcells/GSEA_tcells_ALL.png", plot = last_plot(), device = "png", dpi = 300,
        width = 8, height = 9)

dotplot(go_gsea %>% filter(ONTOLOGY=="BP"), showCategory=30)+ ggtitle("Dotplot for GSEA BP")
ggsave("./figures/GSE144236/tcells/GSEA_tcells_BP.png", plot = last_plot(), device = png, dpi = 300,
       width = 6, height = 9)

dotplot(go_gsea, split = "ONTOLOGY",showCategory=10, x="NES") + 
  facet_grid(ONTOLOGY~., scale="free")+
  ggtitle("Dotplot for GSEA all split")
ggsave("./figures/GSE144236/tcells/GSEA_tcells_ALL_split.png", plot = last_plot(), device = png, dpi = 300,
       width = 8, height = 12)

go_gsea_tbl <- as.tibble(go_gsea)

# KEGG, Reactome, wikipathways

gsea_results_react <- GSEA(
  geneList = ordered_genes, 
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(
    hsa_reactome_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
) 

dotplot(gsea_results_react, x = "NES", showCategory = 30)+ ggtitle("GSEA reactome LFC")
ggsave("./figures/GSE144236/tcells/dotplot_GSEA_reactome_tcells.png", plot = last_plot(), device = png, dpi = 400,
       width = 10, height = 14)

gsea_results_kegg <- GSEA(
  geneList = ordered_genes,
  pvalueCutoff = 0.05, 
  eps = 0,
  seed = TRUE, 
  pAdjustMethod = "BH",
  scoreType = "pos",
  TERM2GENE = dplyr::select(
    hsa_kegg_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_kegg, x = "NES", showCategory = 30)+ ggtitle("GSEA reactome LFC")
ggsave("./figures/GSE144236/tcells/dotplot_GSEA_kegg_tcells.png", plot = last_plot(), device = png, dpi = 400,
       width = 10, height = 8)

# It's still possible to use pathview to extract the significant paths

library("pathview")

keggresids <- c("04060","00190","00480")

foldchanges <- cond.diff.tcell_f$avg_log2FC
names(foldchanges) <- cond.diff.tcell_f$genes
head(foldchanges)
table(is.na(foldchanges))

tmp <- sapply(keggresids, function(pid) pathview(gene.data = foldchanges,
                                                 gene.idtype = "SYMBOL",
                                                 pathway.id = pid,
                                                 species = "hsa",
                                                 kegg.dir="./figures/GSE144236/",
                                                 out.suffix= "_Colored",
                                                 kegg.native = TRUE,
                                                 map.null = FALSE))


gsea_results_wiki <- GSEA(
  geneList = ordered_genes, 
  pvalueCutoff = 0.05, 
  eps = 0, 
  seed = TRUE, 
  pAdjustMethod = "BH", 
  TERM2GENE = dplyr::select(
    hsa_wiki_sets,
    gs_name,
    gene_symbol
  ),
  nPermSimple = 10000
)

dotplot(gsea_results_wiki, x = "NES", showCategory = 30)+ ggtitle("GSEA Wikipathway LFC tcells")

gsea_results_wiki <- as.tibble(gsea_results_wiki)

ggsave("./figures/GSE144236/tcells/dotplot_GSEA_tcells_Wikipathway.png", plot = last_plot(), device = png, dpi = 400)




