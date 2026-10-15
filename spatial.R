library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
plan(multisession(workers = 20))
options(future.globals.maxSize = 8000 * 1024^3)


VlnPlot(xenium.data, features = "nCount_Xenium", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(xenium.data, features = "nCount_Xenium") + theme(legend.position = "right")
wrap_plots(plot1)

ImageDimPlot(xenium.data, fov='fov', molecules = c('SCGB1A1'), nmols = 20000)

ImageFeaturePlot(xenium.data, features= c("ACTA2"), size=0.75, cols=c('white', 'red'))
head(rownames(xenium.data))

cropped.coords <- Crop(xenium.data[["fov"]], x = c(1200, 2900), y = c(3750, 4550), coords = "plot")
xenium.data[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(xenium.data[["zoom"]]) <- "segmentation"
ImageDimPlot(xenium.data, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("ACTA2" , 'CXCL6', 'SFTPC', 'PECAM1'), nmols = 10000)

####Set working Directory
setwd("~/Spatial_Data/CF_gene_markers")

# Find markers for all CF Samples
CF_file <- c('CF19', 'CF20', 'CF27', 'CF34', 'CF35', 'CF40')
load("~/Downloads/Cystic Fibrosis Array/output-XETG00063__0024854__CF19_4-4__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
CF19 <- xenium.data

load("~/Downloads/Cystic Fibrosis Array/output-XETG00063__0024854__CF20_C__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
CF20 <- xenium.data

load("~/Downloads/Cystic Fibrosis Array/output-XETG00063__0024854__CF27_C__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
CF27 <- xenium.data

load("~/Downloads/Cystic Fibrosis Array/output-XETG00063__0024854__CF34_D__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
CF34 <- xenium.data

load("~/Downloads/Cystic Fibrosis Array/output-XETG00063__0024854__CF35_B__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
CF35 <- xenium.data

load("~/Downloads/Cystic Fibrosis Array/output-XETG00063__0024854__CF40_E__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
CF40 <- xenium.data

saveRDS(CF19, 'CF19.rds')
saveRDS(CF20, 'CF20.rds')
saveRDS(CF27, 'CF27.rds')
saveRDS(CF34, 'CF34.rds')
saveRDS(CF35, 'CF35.rds')
saveRDS(CF40, 'CF40.rds')

SpatialDimPlot(CF19)
# Integration 
CF19$lane <- 'CF19'
CF20$lane <- 'C20'
CF27$lane <- 'CF27'
CF34$lane <- 'CF34'
CF35$lane <- 'CF35'
CF40$lane <- 'CF40'

CF_file <- c('CF19', 'CF20', 'CF27', 'CF34', 'CF35', 'CF40')

merged_CF <- merge(CF19, y=c(CF20, CF27, CF34, CF35, CF40),add.cell.ids=CF_file)
VariableFeatures(merged_CF) <- c(VariableFeatures(CF19), VariableFeatures(CF20))
merged_CF <- merged_CF %>%
  RunPCA(npcs = 30, verbose = FALSE)
# one-liner to run Integration
merged_CF = IntegrateLayers(object =merged_CF, method = "RPCAIntegration", normalization.method = "SCT", verbose = T, new.reduction = "integrated.rpca", orig.reduction = "pca")
merged_CF <- merged_CF %>%  
  RunUMAP(reduction = "integrated.rpca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "integrated.rpca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
DimPlot(merged_CF, label=TRUE, split.by = 'lane', raster = FALSE)+ggtitle('CF Samples')
merged_CF_genes <- merged_CF %>% PrepSCTFindMarkers() %>% FindAllMarkers(only.pos = TRUE)
saveRDS(merged_CF, 'merged_CF.csv')
write.csv(merged_CF_genes, 'merged_CF_genes.csv')

################################################################################

# Determine gene markers for all CF samples
CF_list <- c(CF19, CF20, CF27, CF34, CF35, CF40)
for (i in 1:length(CF_list)){
  file <- FindAllMarkers(CF_list[[i]])
  write.csv(file, paste0(CF_file[i], '.csv'))
}

# Find markers for Control Samples
####Set working Directory
setwd("~/Spatial_Data/Ctrl_gene_markers")
Ctrl_file <- c('ID6','MTS7D','MTS12D' ,'MTS29E', 'MTS38D', 'MTS43D')
load("~/Downloads/Non-diseased Array/output-XETG00063__0024799__ID6__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
ID6 <- xenium.data

load("~/Downloads/Non-diseased Array/output-XETG00063__0024799__MTS7D__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
MTS7D <- xenium.data

load("~/Downloads/Non-diseased Array/output-XETG00063__0024799__MTS12D__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
MTS12D <- xenium.data

load("~/Downloads/Non-diseased Array/output-XETG00063__0024799__MTS29E__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
MTS29E <- xenium.data

load("~/Downloads/Non-diseased Array/output-XETG00063__0024799__MTS38D__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
MTS38D <- xenium.data

load("~/Downloads/Non-diseased Array/output-XETG00063__0024799__MTS43D__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
MTS43D <- xenium.data

ID6$lane <- 'ID6'
MTS7D$lane <- 'MTS7D'
MTS12D$lane <- 'MTS12D'
MTS29E$lane <- 'MTS29E'
MTS38D$lane <- 'MTS38D'
MTS43D$lane <- 'MTS43D'

# Determine gene markers for all CF samples
Ctrl_list <- c(ID6,MTS7D, MTS12D, MTS29E, MTS38D, MTS43D)
for (i in 1:length(Ctrl_list)){
  file <- FindAllMarkers(Ctrl_list[[i]], only.pos = TRUE)
  write.csv(file, paste0(Ctrl_file[i], '.csv'))
}

saveRDS(ID6, 'ID6.rds')
saveRDS(MTS12D, 'MTS12D.rds')
saveRDS(MTS29E, 'MTS29E.rds')
saveRDS(MTS38D, 'MTS38D.rds')
saveRDS(MTS43D, 'MTS43D.rds')
saveRDS(MTS7D,'MTS7D')

merged_HC <- merge(ID6, y=c(MTS7D,MTS12D ,MTS29E, MTS38D, MTS43D),add.cell.ids=Ctrl_file)

VariableFeatures(merged_HC) <- c(VariableFeatures(ID6), VariableFeatures(MTS7D), VariableFeatures(MTS12D), VariableFeatures(MTS29E), VariableFeatures(MTS38D), VariableFeatures(MTS43D))
merged_HC <- merged_HC %>%
  RunPCA(npcs = 30, verbose = FALSE)
# one-liner to run Integration
merged_HC = IntegrateLayers(object =merged_HC, method = "RPCAIntegration", normalization.method = "SCT", verbose = T, new.reduction = "integrated.rpca", orig.reduction = "pca")
merged_HC <- merged_HC %>%  
  RunUMAP(reduction = "integrated.rpca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "integrated.rpca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
DimPlot(merged_HC, label=TRUE, split.by = 'lane', raster = FALSE)+ggtitle('HC Samples')
saveRDS(merged_HC, 'merged_HC.rds')

# Find markers for merged_HC object
merged_HC <- readRDS("~/Spatial_Data/Ctrl_gene_markers/merged_HC.rds")
merged_HC_genes <- merged_HC %>% PrepSCTFindMarkers() %>% FindAllMarkers(only.pos = TRUE)
write.csv(merged_HC_genes, 'merged_HC_genes.csv')
# Spatial DimPlot
ImageDimPlot(merged_HC)
ImageFeaturePlot(merged_HC, 'FOXJ1')

###############################################################################
# BOS/Asthma

BOS_file <- c('BOS11A','BOS12D','BOS13E' ,'MTS21A', 'MTS27E', 'MTS43B')
load("~/Downloads/BOS-Asthma array/output-XETG00063__0024851__BOS11A__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
BOS11A <- xenium.data

load("~/Downloads/BOS-Asthma array/output-XETG00063__0024851__BOS12D__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
BOS12D <- xenium.data

load("~/Downloads/BOS-Asthma array/output-XETG00063__0024851__BOS13E__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
BOS13E <- xenium.data

load("~/Downloads/BOS-Asthma array/output-XETG00063__0024851__MTS21A__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
MTS21A <- xenium.data

load("~/Downloads/BOS-Asthma array/output-XETG00063__0024851__MTS27E__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
MTS27E <- xenium.data

load("~/Downloads/BOS-Asthma array/output-XETG00063__0024851__MTS48B__20240711__193926/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
MTS43B <- xenium.data

BOS_list <- c(BOS11A, BOS12D, BOS13E, MTS21A, MTS27E, MTS43B)
for (i in 1:length(BOS_list)){
  file <- FindAllMarkers(BOS_list[[i]])
  write.csv(file, paste0(BOS_file[i], '.csv'))
}
saveRDS(BOS11A, 'BOS11A.rds')
saveRDS(BOS12D, 'BOS12D.rds')
saveRDS(BOS13E, 'BOS13E.rds')
saveRDS(MTS21A, 'MTS21A.rds')
saveRDS(MTS27E, 'MTS27E.rds')
saveRDS(MTS43B, 'MTS43B.rds')

################################################################################
# COPD

COPD_file <- c('C48A','C50C','C51E' ,'C52E', 'C53E', 'C56E')
load("~/Downloads/COPD array/output-XETG00063__0024796__C48A__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
C48A <- xenium.data

load("~/Downloads/COPD array/output-XETG00063__0024796__C50C__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
C50C <- xenium.data

load("~/Downloads/COPD array/output-XETG00063__0024796__C51E__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
C51E <- xenium.data

load("~/Downloads/COPD array/output-XETG00063__0024796__C52E__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
C52E <- xenium.data

load("~/Downloads/COPD array/output-XETG00063__0024796__C53E__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
C53E <- xenium.data

load("~/Downloads/COPD array/output-XETG00063__0024796__C56E__20240703__220049/Seurat_Cluster_Analysis_0.3_Resolution/Seurat.RData")
C56E <- xenium.data

COPD_list <- c(C48A,C50C,C51E ,C52E, C53E, C56E)
for (i in 1:length(COPD_list)){
  file <- FindAllMarkers(COPD_list[[i]])
  write.csv(file, paste0(COPD_file[i], '.csv'))
}
saveRDS(C48A, 'C48A.rds')
saveRDS(C50C, 'C50C.rds')
saveRDS(C51E, 'C51E.rds')
saveRDS(C52E, 'C52E.rds')
saveRDS(C53E, 'C53E.rds')
saveRDS(C56E, 'C56E.rds')
