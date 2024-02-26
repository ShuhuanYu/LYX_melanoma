library(DoubletFinder)
library(Seurat)
library(hdf5r)
library(ggplot2)

samples=c("CR177","CR178","CR186")

### pre-processing data  
removeDoublet.cellNumber <- matrix(rep(0,3),3,3)
rownames(removeDoublet.cellNumber) <- samples
colnames(removeDoublet.cellNumber) <- c("raw_cell_numbers","cell_numbers_after_qc","cell_numbers_after_removedoublets")

for(i in 2:length(samples)){  
  sample <- samples[i]
  print(paste0("Now processing: ",sample))

  h5_file <- paste0("./data/",sample,"/filtered_feature_bc_matrix.h5")
  inputdata.10x <- Read10X_h5(h5_file)
  mouse_melonoma <- CreateSeuratObject(counts = inputdata.10x, project = "melonoma", min.cells = 3, min.features = 200)
  mouse_melonoma[["percent.mt"]] <- PercentageFeatureSet(mouse_melonoma, pattern = "^MT-")
  
  raw_cell_number <- dim(mouse_melonoma)[2]

  #quality control
  mouse_melonoma <- subset(mouse_melonoma, subset = nFeature_RNA > 200 & nFeature_RNA < 25000 & percent.mt < 5)
  
  cell_numbers_after_qc <- dim(mouse_melonoma)[2]

  mouse_melonoma <- NormalizeData(mouse_melonoma)
  mouse_melonoma <- FindVariableFeatures(mouse_melonoma, selection.method = "vst", nfeatures = 2000)
  mouse_melonoma <- ScaleData(mouse_melonoma)
  mouse_melonoma <- RunPCA(mouse_melonoma)
  mouse_melonoma <- RunUMAP(mouse_melonoma, dims = 1:20)

  # cluster cells
  mouse_melonoma <- FindNeighbors(mouse_melonoma, dims = 1:20)
  mouse_melonoma <- FindClusters(mouse_melonoma,resolution=0.2)
  
  # pK Identification (no ground-truth) 
  sweep.res.list <- paramSweep(mouse_melonoma, PCs = 1:20, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK<-as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  # remofve doublets by DoubletFinder
  mouse_melonoma$celltype <- Idents(mouse_melonoma)
  DoubletRate=ncol(mouse_melonoma)*8*1e-6
  homotypic.prop=modelHomotypic(mouse_melonoma@meta.data$celltype)
  nExp_poi=round(DoubletRate*length(mouse_melonoma$celltype))
  nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
  mouse_melonoma=doubletFinder(mouse_melonoma, PCs = 1:20, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = FALSE)
  DF.clarrsications=names(mouse_melonoma@meta.data)[grep("DF.classifications",names(mouse_melonoma@meta.data))]
  
  # Plot results
  p1 <- DimPlot(mouse_melonoma, reduction = "umap", group.by = DF.clarrsications) + ggtitle("Doublets")
  ggsave(paste0("./plots/",sample,"_DoubletFinder.pdf"),p1,width = 10,height = 10)

  # filter doublets
  Idents(mouse_melonoma ) <- DF.clarrsications
  mouse_melonoma<-subset(x = mouse_melonoma, idents="Singlet")   #过滤细胞，只保留经过检测的Singlet
  
  filtered_cell_number <- dim(mouse_melonoma)[2]

  removeDoublet.cellNumber[i,1] = raw_cell_number
  removeDoublet.cellNumber[i,2] = cell_numbers_after_qc
  removeDoublet.cellNumber[i,3] = filtered_cell_number

  print(paste0("raw cell numbers: ",raw_cell_number))
  print(paste0("cell numbers after quality control: ", cell_numbers_after_qc))
  print(paste0("filtered cell numbers: ",filtered_cell_number))

  # save filtered rna_counts 
  filtered_barcodes <- colnames(mouse_melonoma[["RNA"]]$counts)
  filtered.rna.counts = mouse_melonoma[["RNA"]]$counts[,filtered_barcodes]
  colnames(filtered.rna.counts) = paste0(sample,".",colnames(filtered.rna.counts))

  write.table(filtered.rna.counts, paste0("./tables/","/",sample,"_filtered.rna.counts.txt"),row.names = T,sep = "\t",quote = F)

  # rm(list = ls())
  gc()
}

write.table(removeDoublet.cellNumber,"./tables/removeDoublet.details.txt",sep="\t",quote=F,row.names = T)

# Session info:
# R version 4.3.1 (2023-06-16 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 11 x64 (build 22621)

# Matrix products: default


# locale:
# [1] LC_COLLATE=Chinese (Simplified)_China.utf8 
# [2] LC_CTYPE=Chinese (Simplified)_China.utf8
# [3] LC_MONETARY=Chinese (Simplified)_China.utf8
# [4] LC_NUMERIC=C
# [5] LC_TIME=Chinese (Simplified)_China.utf8

# time zone: Asia/Shanghai
# tzcode source: internal

# attached base packages:
# [1] parallel  stats     graphics  grDevices utils     datasets  methods
# [8] base

# other attached packages:
#  [1] ggplot2_3.4.4       ROCR_1.0-11         KernSmooth_2.23-22
#  [4] fields_15.2         viridisLite_0.4.2   spam_2.10-0
#  [7] hdf5r_1.3.9         Seurat_5.0.1        SeuratObject_5.0.1
# [10] sp_2.1-2            DoubletFinder_2.0.4

# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3     jsonlite_1.8.8         magrittr_2.0.3
#   [4] spatstat.utils_3.0-4   farver_2.1.1           ragg_1.2.7
#   [7] vctrs_0.6.5            spatstat.explore_3.2-5 htmltools_0.5.7       
#  [10] curl_5.2.0             sctransform_0.4.1      parallelly_1.36.0
#  [13] htmlwidgets_1.6.4      desc_1.4.3             ica_1.0-3
#  [16] plyr_1.8.9             plotly_4.10.3          zoo_1.8-12
#  [19] igraph_1.6.0           mime_0.12              lifecycle_1.0.4
#  [22] pkgconfig_2.0.3        Matrix_1.6-4           R6_2.5.1
#  [25] fastmap_1.1.1          fitdistrplus_1.1-11    future_1.33.1
#  [28] shiny_1.8.0            digest_0.6.33          colorspace_2.1-0
#  [31] patchwork_1.1.3        ps_1.7.5               tensor_1.5
#  [34] RSpectra_0.16-1        irlba_2.3.5.1          textshaping_0.3.7     
#  [37] labeling_0.4.3         progressr_0.14.0       fansi_1.0.6
#  [40] spatstat.sparse_3.0-3  httr_1.4.7             polyclip_1.10-6
#  [43] abind_1.4-5            compiler_4.3.1         remotes_2.4.2.1
#  [46] bit64_4.0.5            withr_2.5.2            fastDummies_1.7.3
#  [49] pkgbuild_1.4.3         maps_3.4.2             MASS_7.3-60
#  [52] tools_4.3.1            lmtest_0.9-40          httpuv_1.6.13
#  [55] future.apply_1.11.1    goftest_1.2-3          glue_1.6.2
#  [58] callr_3.7.3            nlme_3.1-164           promises_1.2.1
#  [61] grid_4.3.1             Rtsne_0.17             cluster_2.1.6
#  [64] reshape2_1.4.4         generics_0.1.3         gtable_0.3.4
#  [67] spatstat.data_3.0-3    tidyr_1.3.0            data.table_1.14.10
#  [70] utf8_1.2.4             spatstat.geom_3.2-7    RcppAnnoy_0.0.21
#  [73] ggrepel_0.9.4          RANN_2.6.1             pillar_1.9.0
#  [76] stringr_1.5.1          RcppHNSW_0.5.0         later_1.3.2
#  [79] splines_4.3.1          dplyr_1.1.4            lattice_0.22-5
#  [82] survival_3.5-7         bit_4.0.5              deldir_2.0-2
#  [85] tidyselect_1.2.0       miniUI_0.1.1.1         pbapply_1.7-2
#  [88] gridExtra_2.3          scattermore_1.2        matrixStats_1.2.0
#  [91] stringi_1.8.3          lazyeval_0.2.2         codetools_0.2-19
#  [94] tibble_3.2.1           cli_3.6.2              uwot_0.1.16
#  [97] xtable_1.8-4           reticulate_1.34.0      systemfonts_1.0.5
# [100] munsell_0.5.0          processx_3.8.3         Rcpp_1.0.12
# [103] globals_0.16.2         spatstat.random_3.2-2  png_0.1-8
# [106] ellipsis_0.3.2         dotCall64_1.1-1        listenv_0.9.0
# [109] scales_1.3.0           ggridges_0.5.5         leiden_0.4.3.1
# [112] purrr_1.0.2            rlang_1.1.2            cowplot_1.1.2