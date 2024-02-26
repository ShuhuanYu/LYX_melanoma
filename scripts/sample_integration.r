library(Seurat)
library(EnsDb.Mmusculus.v79)
library(Signac)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(patchwork)
library(hdf5r)
library(harmony)
library(cowplot)
library(SeuratWrappers)
library(scCustomize)

### Integration BD samples

CR177.counts <- read.table("./tables/CR177_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
CR178.counts <- read.table("./tables/CR178_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)
CR186.counts <- read.table("./tables/CR186_filtered.rna.counts.txt", header = T,sep = "\t",row.names = 1)

genes <- intersect(intersect(rownames(CR177.counts), rownames(CR178.counts)), rownames(CR186.counts))

CR177.counts <- CR177.counts[genes, ]
CR178.counts <- CR178.counts[genes, ]
CR186.counts <- CR186.counts[genes, ]

all.samples.counts <- cbind(CR177.counts,CR178.counts,CR186.counts)
rm("CR177.counts", "CR178.counts", "CR186.counts")

melanoma <- CreateSeuratObject(counts = all.samples.counts, project = "melanoma", min.cells = 0, min.features = 0) %>%
    Seurat::NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
    ScaleData(verbose = FALSE) %>% 
    RunPCA(pc.genes = melanoma@var.genes, npcs = 20, verbose = FALSE)

melanoma@meta.data$Samples <- c(rep("CR177", 5908), rep("CR178", 6741), rep("CR186", 8044))

# run Harmony
options(repr.plot.height = 2.5, repr.plot.width = 6)
melanoma <- melanoma %>% 
    RunHarmony("Samples", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(melanoma, 'harmony')
harmony_embeddings[1:5, 1:5]

options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = melanoma, reduction = "harmony", pt.size = .1, group.by = "Samples")
p2 <- VlnPlot(object = melanoma, features = "harmony_1", group.by = "Samples", pt.size = .1)
ggsave("./plots/harmony_results.pdf", p1 + p2, width = 15, height = 8)


# UMAP analysis
melanoma <- RunUMAP(melanoma, reduction = "harmony", dims = 1:20)
melanoma <- FindNeighbors(melanoma, reduction = "harmony", dims = 1:20)
melanoma <- FindClusters(melanoma, resolution = 1.5) 
# ElbowPlot(melanoma, ndims = 20)
table(Idents(melanoma))

p1 <- DimPlot(melanoma, reduction = "umap", group.by = "Samples", pt.size = .1, split.by = 'Samples')
ggsave("./plots/plot_unlabel_Harmony_dimplot.pdf", p1, width = 20,height = 7)

p2 <- DimPlot(melanoma, reduction = "umap", label = TRUE, pt.size = .1)
ggsave("./plots/plot_unlabel_Harmony_umap.pdf", p2, width = 7,height = 7)

# annotate manually
marker.genes <- c("Cd79a", "Ms4a1", "Cd79b", "Cd19", # B_cells
                    "Cd74", "Irf8", "Itgax", "Clec9a", "H2-Ab1", # Dendric_cells
                    "Vwf", "Pecam1", "Cldn5", "Flt1", "Cdh5", # Endothelial
                    "Epcam", "Krt5", "Krt14", "Krt19", # Epithelial
                    "Dcn", "Col1a1", "Pdgfra", # Fibroblast
                    # "Ptprc", # Immune_cells
                    "Aif1", "Csf1r", "C1qb", "C1qa", # Macrophage
                    # "Kit", # Mast_cell
                    "Mlana", # Melanocyte
                    "S100a8", "S100a9", "Csf3r", # Neutrophil
                    "Nkg7", "Klrd1", "Gzma", "Prf1", # NK_cell
                    "Cd79a", "Jchain", "Tnfrsf17", # Plasma
                    "Rpe65", # RPE
                    # "Plp1", "Ncmap", "Scn7a", # Schwann
                    # "Myh11", "Acta2", "Myl9", "Tagln", # Smooth_muscle_cell
                    "Cd2", "Cd3d", "Cd3e", "Cd3g", "Trac" # T_cell
                 )
cd_genes=unique(marker.genes)
#
cell.type.cols <- c("#DC143C","#20B2AA","#FFA500","#9370DB","#F08080","#1E90FF",
            "#808000","#CCCCFF","#f04f04","#800080","#D2691E",
            "#87CEEB","#40E0D0","#008B8B","#228B22","#E9967A",
            "#4682B4","#F0E68C","#EE82EE","#FF6347","#8B4513",
            "#35c041","#b1ef45","#cff819","#0c0fae","#530f31","#f41f8a",
            "#1f4723","#91938e","#364202","#06074c","#280215","#37031d",
            "#7ef688"
            )
p1 <- Stacked_VlnPlot(melanoma,features = cd_genes,colors_use = cell.type.cols)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
p2 <- DotPlot(melanoma, features = cd_genes)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("./plots/seurat.anno_marker.genes.plot.pdf", p2+p1, width = 30, height = 25)
ggsave("./plots/seurat.anno_marker.genes.dotplot.pdf", p2, width = 20, height = 10)

melanoma <- RenameIdents(object = melanoma, 
                    "0" = "Melanocyte or Melanoma", 
                    "1" = "Melanocyte or Melanoma",
                    "2" = "T cell",
                    "3" = "Macrophage",
                    "4" = "Neutrophil",
                    "5" = "Melanocyte or Melanoma",
                    "6" = "Macrophage", 
                    "7" = "Melanocyte or Melanoma", 
                    "8" = "Melanocyte or Melanoma",
                    "9" = "T cell", 
                    "10" = "Melanocyte or Melanoma", 
                    "11" = "Melanocyte or Melanoma", 
                    "12" = "Melanocyte or Melanoma",
                    "13" = "Melanocyte or Melanoma", 
                    "14" = "Melanocyte or Melanoma",
                    "15" = "Macrophage",
                    "16" = "T cell",
                    "17" = "Macrophage", 
                    "18" = "T cell", 
                    "19" = "Melanocyte or Melanoma",
                    "20" = "Melanocyte or Melanoma",
                    "21" = "DC",
                    "22" = "Fibroblast",
                    "23" = "Melanocyte or Melanoma",
                    "24" = "Melanocyte or Melanoma",
                    "25" = "RPE",
                    '26' = 'NK cell',
                    '27' = 'B cell',
                    '28' = 'DC',
                    '29' = 'Endothelial',
                    '30' = 'Melanocyte or Melanoma',
                    '31' = 'Melanocyte or Melanoma',
                    '32' = 'Melanocyte or Melanoma',
                    '33' = 'Epithelial'
                    )

melanoma$celltype <- Idents(melanoma)

p1 <- DimPlot(melanoma, reduction = "umap",label=TRUE, cols = cell.type.cols[1:11])
ggsave("./plots/Plots-label-umap.pdf",p1, width = 10, height = 8)

levels(melanoma) <- c("B cell", "DC", "Endothelial", "Epithelial",
                      "Fibroblast", "Macrophage", "Melanocyte or Melanoma", "Neutrophil",
                      "NK cell", "RPE", "T cell"
                      )
p2 <- DotPlot(melanoma, features = rev(cd_genes)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave("./plots/Plots-label-dotplot.pdf",p2, width = 15, height = 8)


saveRDS(melanoma,"./tables/melanoma_scRNA_harmony.rds")
# melanoma <- readRDS("./tables/melanoma_scRNA_harmony.rds")
