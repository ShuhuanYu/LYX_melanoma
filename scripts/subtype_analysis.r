library(Seurat)
library(ggplot2)
library(gplots)
library(patchwork)
library(dplyr)
library(clusterProfiler)
library(org.Mmu.eg.db)
library(harmony)
library(cowplot)
library(scCustomize)

# load data
melanoma <- readRDS("./tables/melanoma_scRNA_harmony.rds")

# extract celltype-specific cells
celltype <- "T cell" # 3206 cells

celltype_cells <- melanoma[, Idents(melanoma) %in% celltype] 

celltype_cells <- NormalizeData(celltype_cells, normalization.method = "LogNormalize", scale.factor = 1e4) 
celltype_cells <- FindVariableFeatures(celltype_cells, selection.method = 'vst', nfeatures = 2000)
celltype_cells <- ScaleData(celltype_cells,features=rownames(celltype_cells))
celltype_cells <- RunPCA(celltype_cells, features = VariableFeatures(object = celltype_cells)) 
ElbowPlot(celltype_cells, ndims=20, reduction="pca")

celltype_cells <- RunHarmony(celltype_cells, group.by.vars="Samples", plot_convergence=TRUE)

celltype_cells <- FindNeighbors(celltype_cells, dims = 1:20)
celltype_cells <- FindClusters(celltype_cells, resolution = 2.0 )
table(celltype_cells$seurat_clusters) 

celltype_cells <- RunUMAP(celltype_cells, reduction = "harmony", dims = 1:20)

p1 <- DimPlot(celltype_cells, reduction = "umap", label = TRUE, pt.size = .1)
ggsave("./plots/T_cell_unlabel_umap.pdf", p1, width = 10,height = 8)

# cell annotation
rownames(celltype_cells) <- toupper(rownames(celltype_cells))
marker_genes <- c(
    "FOXP3", "TNFRSF4", "IKZF2", "IL2RA", # T regulatory cell
    "MAF", "CD200", "GNG4", "CHN1", "IGFL2", "ITM2A", "CPM", "NR3C1", # Follicular T cell
    "IL7R", "TCF7", "SELL", "LEF1", "CCR7", # Naive T cell
    "CD8A", "ZNF683", # CD8+ T effector memory cell
    "KLRK1", "ITGAE", "CD103", # CD8+ resident memory cell
    "I17R", "CD4", # CD4+ T cell
    "TRDC", "TRDG2", "TRDC1", "TRDC2", # CD8+ gamma delta T cell
    # "MKI67", "STMN1", "HMGB2", "TUBB", "TUBA1B", # Mitotic CD8+ T cell
    "PRF1", "GZMA", "GZMK", "NKG7", # Cytotoxic CD8+ T cell
    # "IRF4", "CREM", "NR4A2", # T helper 17
    "STAT4", "IFNG", "IL12RB2" # T helper 1
    # "GATA3", "STAT6", "IL4", # T helper 2
    # "MAF", "CXCR5", "PDCD1", "CXCL13", # T follocular helper
    # "XCL1", "FCGR3A", "KLRD1", "KLRF1" # NK cell
    )
#     "IL2RA", "FOXP3", "IKZF2", "TGFB1", "TGFB3", "TGFBI", "TGFBR1", # Treg

# )
cd_genes <- unique(marker_genes)

p2 <- DotPlot(celltype_cells, features = cd_genes)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("./plots/T_cell_unlabel_dotplot.pdf", p2, width = 12,height = 7)

xxxmarkers<- FindMarkers(celltype_cells, ident.1 = "4", only.pos = TRUE)
View(xxxmarkers)

celltype_cells <- RenameIdents(celltype_cells, 
                        '0' = 'Follicular T cell',
                        '1' = 'Cytotoxic CD8+ T cell',
                        '2' = 'Cytotoxic CD8+ T cell',
                        '3' = 'Cytotoxic CD8+ T cell',
                        '4' = 'Naive T cell',
                        '5' = 'Cytotoxic CD8+ T cell',
                        '6' = 'Follicular T cell',
                        '7' = 'Cytotoxic CD8+ T cell',
                        '8' = 'T helper 1',
                        '9' = 'Cytotoxic CD8+ T cell',
                        '10' = 'Follicular T cell',
                        '11' = 'Follicular T cell',
                        '12' = 'Cytotoxic CD8+ T cell',
                        '13' = 'T regulatory cell',
                        '14' = 'Naive T cell',
                        '15' = 'Cytotoxic CD8+ T cell',
                        '16' = 'Follicular T cell'
                        )
celltype_cells$celltype <- Idents(celltype_cells)

p1 <- DimPlot(celltype_cells, reduction ="umap", label = TRUE)
p2 <- DimPlot(celltype_cells, reduction="umap", group.by="Samples")
ggsave("./plots/T_cell_subcluster_umap.pdf", p1,width=8,height=6)

saveRDS(celltype_cells,"./tables/T_cell_subclusters.rds")

### depletion factor ###
dep_factors <- c("LAG3", "PD1", "CTLA4", "TIGIT", "HAVCR2", "TIM3", "TNFRSF9", "4-1BB")

p = FeaturePlot(celltype_cells, features = dep_factors, reduction = "umap")
ggsave("./plots/CD8+T_cell_depletion_factors_featureplot.pdf", p, width = 10, height = 10)
