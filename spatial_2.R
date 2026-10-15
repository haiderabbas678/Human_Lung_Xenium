library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
plan(multisession(workers = 20))
options(future.globals.maxSize = 8000 * 1024^3)
setwd("~/Spatial_Data")
################# CF Samples QC and pre-processing##############################
CF19 <- readRDS("~/Spatial_Data/CF/CF19.rds")
CF20 <- readRDS("~/Spatial_Data/CF/CF20.rds")
CF27 <- readRDS("~/Spatial_Data/CF/CF27.rds")
CF35 <- readRDS("~/Spatial_Data/CF/CF35.rds")
CF40 <- readRDS("~/Spatial_Data/CF/CF40.rds")

# Removing low quality labelled clusters
CF19 <- subset(CF19, idents = 'LQ', invert=TRUE)
CF19$Group <- 'CF'
CF19$ID <- 'CF19'

CF20 <- subset(CF20, idents = 'LQ', invert=TRUE)
CF20$Group <- 'CF'
CF20$ID <- 'CF20'

CF27$Group <- 'CF'
CF27$ID <- 'CF27'

CF35 <- subset(CF35, idents = 'LQ', invert=TRUE)
CF35$Group <- 'CF'
CF35$ID <- 'CF35'

CF40 <- subset(CF40, idents = c('LQ', '?LQ'), invert=TRUE)
CF40$Group <- 'CF'
CF40$ID <- 'CF40'

# nFeature_Xenium Violin Plots
#CF_list <- c(CF19, CF20, CF27, CF35, CF40)
#CF_file <- c('CF19', 'CF20', 'CF27', 'CF35', 'CF40')
#for (i in 1:length(CF_list)){
 # violin <- VlnPlot(CF_list[[i]], features = c("nFeature_Xenium"), raster = FALSE, pt.size = 0.001)
  #pdf(paste0(CF_file[i], '.pdf'),8,8)
  #plot(violin)
  #dev.off()
#}
################# BOS Samples QC and pre-processing#############################
BOS11A <- readRDS("~/Spatial_Data/BOS/BOS11A_sub.rds")
BOS12D <- readRDS("~/Spatial_Data/BOS/BOS12D_sub.rds")
BOS13E <- readRDS("~/Spatial_Data/BOS/BOS13E_sub.rds")
MTS21A <- readRDS("~/Spatial_Data/BOS/MTS21A_sub.rds")
MTS48B <- readRDS("~/Spatial_Data/BOS/MTS48B_sub.rds")
final_obj$ID[final_obj$ID %in% c('MTS43B')] <- paste0('MTS48B')
levels(factor(final_obj$ID))
MTS27E <- readRDS("~/Spatial_Data/BOS/MTS27E.rds")

BOS11A$Group <- 'BOS'
BOS11A$ID <- 'BOS11A'

BOS12D$Group <- 'BOS'
BOS12D$ID <- 'BOS12D'

BOS13E$Group <- 'BOS'
BOS13E$ID <- 'BOS13E'

MTS21A$Group <- 'Asthma'
MTS21A$ID <- 'MTS21A'

MTS43B$Group <- 'Asthma'
MTS48B$ID <- 'MTS48B'

MTS27E$Group <- 'Asthma'
MTS27E$ID <- 'MTS27E'
# nFeature_Xenium Violin Plots
#BOS_list <- c(BOS11A, BOS12D, BOS13E, MTS21A, MTS43B, MTS27E)
#BOS_file <- c('BOS11A', 'BOS12D', 'BOS13E', 'MTS21A', 'MTS43B', 'MTS27E')
#for (i in 1:length(BOS_list)){
  #violin <- VlnPlot(BOS_list[[i]], features = c("nFeature_Xenium"), raster = FALSE, pt.size = 0.001)
  #pdf(paste0(BOS_file[i], '.pdf'),8,8)
  #plot(violin)
  #dev.off()
#}
################# HC Samples QC and pre-processing##############################
ID6 <- readRDS("~/Spatial_Data/HC/ID6.rds")
MTS7D <- readRDS("~/Spatial_Data/HC/MTS7D.rds")
MTS12D <- readRDS("~/Spatial_Data/HC/MTS12D.rds")
MTS29E <- readRDS("~/Spatial_Data/HC/MTS29E.rds")
MTS38D <- readRDS("~/Spatial_Data/HC/MTS38D.rds")
MTS43D <- readRDS("~/Spatial_Data/HC/MTS43D.rds")

# Removing low quality labelled clusters
ID6 <- subset(ID6, idents = 'LQ', invert=TRUE)
ID6$Group <- 'HC'
ID6$ID <- 'ID6'

MTS7D <- subset(MTS7D, idents = c('LQ_PMN', 'LQ_FB'), invert=TRUE)
MTS7D$Group <- 'HC'
MTS7D$ID <- 'MTS7D'

MTS12D <- subset(MTS12D, idents = 'LQ', invert=TRUE)
MTS12D$Group <- 'HC'
MTS12D$ID <- 'MTS12D'

MTS29E$Group <- 'HC'
MTS29E$ID <- 'MTS29E'

MTS38D <- subset(MTS38D, idents = 'LQ_SMC', invert=TRUE)
MTS38D$Group <- 'HC'
MTS38D$ID <- 'MTS38D'

MTS43D$Group <- 'HC'
MTS43D$ID <- 'MTS43D'
# nFeature_Xenium Violin Plots
#Ctrl_list <- c(ID6,MTS7D, MTS12D, MTS29E, MTS38D, MTS43D)
#Ctrl_file <- c('ID6','MTS7D', 'MTS12D', 'MTS29E', 'MTS38D', 'MTS43D')
#for (i in 1:length(Ctrl_list)){
 # violin <- VlnPlot(Ctrl_list[[i]], features = c("nFeature_Xenium"), raster = FALSE, pt.size = 0.001)
  #pdf(paste0(Ctrl_file[i], '.pdf'),8,8)
  #plot(violin)
  #dev.off()
#}



################# COPD Samples QC and pre-processing############################
C48A <- readRDS("~/Spatial_Data/COPD/C48A.rds")
C50C <- readRDS("~/Spatial_Data/COPD/C50C.rds")
C51E <- readRDS("~/Spatial_Data/COPD/C51E.rds")
C52E <- readRDS("~/Spatial_Data/COPD/C52E.rds")
C56E <- readRDS("~/Spatial_Data/COPD/C56E.rds")

# Removing low quality labelled clusters
C48A <- subset(C48A, idents = 'LQ', invert=TRUE)
C48A$Group <- 'COPD'
C48A$ID <- 'C48A'

C50C$Group <- 'COPD'
C50C$ID <- 'C50C'

C51E <- subset(C51E, idents = 'LQ', invert=TRUE)
C51E$Group <- 'COPD'
C51E$ID <- 'C51E'

C52E <- subset(C52E, idents = 'LQ', invert=TRUE)
C52E$Group <- 'COPD'
C52E$ID <- 'C52E'

C56E$Group <- 'COPD'
C56E$ID <- 'C56E'

#COPD_list <- c(C48A,C50C,C51E ,C52E, C56E)
#COPD_file <- c('C48A','C50C','C51E' ,'C52E', 'C56E')
#for (i in 1:length(COPD_list)){
 # violin <- VlnPlot(COPD_list[[i]], features = c("nFeature_Xenium"), raster = FALSE, pt.size=0.001)
  #pdf(paste0(COPD_file[i], '.pdf'),8,8)
  #plot(violin)
  #dev.off()
#}

#VlnPlot(C48A, features = c("nFeature_Xenium"),pt.size=0.001, raster = FALSE)

################# Integration of Data Sets######################################
cell.ids <- c('ID6','MTS7D', 'MTS12D', 'MTS29E', 'MTS38D', 'MTS43D','C48A','C50C','C51E' ,'C52E', 'C56E')
merged_COPD_HC <- merge(ID6, y=c(MTS7D, MTS12D, MTS29E, MTS38D, MTS43D, C48A,C50C,C51E ,C52E, C56E),add.cell.ids=cell.ids)
VariableFeatures(merged_COPD_HC) <- c(VariableFeatures(ID6), VariableFeatures(MTS7D), VariableFeatures(MTS12D), VariableFeatures(MTS29E), VariableFeatures(MTS38D), VariableFeatures(MTS43D), VariableFeatures(C48A),
                                      VariableFeatures(C50C), VariableFeatures(C51E), VariableFeatures(C52E), VariableFeatures(C56E))
merged_COPD_HC <- merged_COPD_HC %>%
  RunPCA(npcs = 30, verbose = FALSE)
# one-liner to run Integration
merged_COPD_HC = IntegrateLayers(object =merged_COPD_HC, method = "HarmonyIntegration", normalization.method = "SCT", verbose = T, new.reduction = "harmony", orig.reduction = "pca")
merged_COPD_HC <- merged_COPD_HC %>%  
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

saveRDS(merged_COPD_HC, 'merged_COPD_HC.rds')
merged_COPD_HC_genes <- merged_COPD_HC %>% PrepSCTFindMarkers() %>% FindAllMarkers(only.pos = TRUE)
write.csv(merged_COPD_HC_genes, 'merged_COPD_HC_genes.csv')
################################################################################
cell.ids <- c('ID6','MTS7D', 'MTS12D', 'MTS29E', 'MTS38D', 'MTS43D','C48A','C50C','C51E' ,'C52E', 'C56E', 'CF19', 'CF20', 'CF27', 'CF35', 'CF40')
merged_final <- merge(ID6, y=c(MTS7D, MTS12D, MTS29E, MTS38D, MTS43D, C48A,C50C,C51E ,C52E, C56E, CF19, CF20, CF27, CF35, CF40),add.cell.ids=cell.ids)
VariableFeatures(merged_final) <- c(VariableFeatures(ID6), VariableFeatures(MTS7D), VariableFeatures(MTS12D), VariableFeatures(MTS29E), VariableFeatures(MTS38D), VariableFeatures(MTS43D), VariableFeatures(C48A),
                                      VariableFeatures(C50C), VariableFeatures(C51E), VariableFeatures(C52E), VariableFeatures(C56E), VariableFeatures(CF19), VariableFeatures(CF20), VariableFeatures(CF27), VariableFeatures(CF35),
                                    VariableFeatures(CF40))
merged_final <- merged_final %>%
  RunPCA(npcs = 30, verbose = FALSE)
# one-liner to run Integration
merged_final = IntegrateLayers(object =merged_final, method = "HarmonyIntegration", normalization.method = "SCT", verbose = T, new.reduction = "harmony", orig.reduction = "pca")
merged_final <- merged_final %>%  
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

DimPlot(merged_final, label = TRUE, raster=FALSE)

saveRDS(merged_final, 'merged_final.rds')
merged_final_genes <- merged_final %>% PrepSCTFindMarkers() %>% FindAllMarkers(only.pos = TRUE)
write.csv(merged_final_genes, 'merged_final_genes.csv')

FeaturePlot(merged_final,c('SFTPD'),  cols = c('grey94', 'darkblue'), order=T)s
#Table (cells vs condition)
final_table <- table(final_obj$Cell_Class, final_obj$Group)
final_table[,1] <- final_table[,1]/sum(final_table[,1])*100
final_table[,2] <- final_table[,2]/sum(final_table[,2])*100
final_table[,3] <- final_table[,3]/sum(final_table[,3])*100
final_table

################# Combined Object Inetgration ##################################
cell.ids <- c('ID6','MTS7D', 'MTS12D', 'MTS29E', 'MTS38D', 'MTS43D','C48A','C50C','C51E' ,'C52E', 'C56E', 'CF19', 'CF20', 'CF27', 'CF35', 'CF40', 'BOS11A', 'BOS12D', 'BOS13E', 'MTS21A', 'MTS48B', 'MTS27E')
comb_final <- merge(ID6, y=c(MTS7D, MTS12D, MTS29E, MTS38D, MTS43D, C48A,C50C,C51E ,C52E, C56E, CF19, CF20, CF27, CF35, CF40, BOS11A, BOS12D, BOS13E, MTS21A, MTS43B, MTS27E),add.cell.ids=cell.ids)
VariableFeatures(comb_final) <- c(VariableFeatures(ID6), VariableFeatures(MTS7D), VariableFeatures(MTS12D), VariableFeatures(MTS29E), VariableFeatures(MTS38D), VariableFeatures(MTS43D), VariableFeatures(C48A),
                                    VariableFeatures(C50C), VariableFeatures(C51E), VariableFeatures(C52E), VariableFeatures(C56E), VariableFeatures(CF19), VariableFeatures(CF20), VariableFeatures(CF27), VariableFeatures(CF35),
                                    VariableFeatures(CF40),  VariableFeatures(BOS11A),  VariableFeatures(BOS12D),  VariableFeatures(BOS13E),  VariableFeatures(MTS21A),  VariableFeatures(MTS43B), VariableFeatures(MTS27E))
comb_final <- comb_final %>%
  RunPCA(npcs = 30, verbose = FALSE)
ElbowPlot(comb_final)
# one-liner to run Integration
comb_final = IntegrateLayers(object =comb_final, method = "HarmonyIntegration", normalization.method = "SCT", verbose = T, new.reduction = "harmony", orig.reduction = "pca")
comb_final <- comb_final %>%  
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
saveRDS(comb_final, 'comb_final.rds')


comb_final_genes <- comb_final %>% PrepSCTFindMarkers() %>% FindAllMarkers(only.pos = TRUE)
saveRDS(comb_final, 'comb_final.rds')
write.csv(comb_final_genes, 'comb_final_genes.csv')
DimPlot(comb_final, label=TRUE)
VlnPlot(comb_final, features = 'nFeature_Xenium', pt.size = 0)
FeaturePlot(comb_final, features = c("MMP12"), cols = c('grey94', 'darkblue'))

# Annotation
cell.idents <- c('AT2', 'Mono/Macs', 'gCap', 'FCN3+_EC', 'Cil', 'SMC', 'FB', 'AT1',
                 'ThC', 'Basal', 'APOD+_FB', 'GB', 'TkC', 'PLA2G2A+_FB', 'AM', 'DES+_SMC', 
                 'aCap', 'Perivasc_macs', 'Gob', 'Peri', 'PMN', 'Sec', 'BC', 'AT2_EC_DB', 
                 'Mast', 'LEC', 'DMBT1+_AT2', 'CXCL+_Immune', 'LQ', 'cDC', 'LQ', 'LQ', 'Basal', 
                 'AT2', 'LQ', 'LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ',
                 'LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ','LQ')
comb_final@active.ident <- comb_final$seurat_clusters
names(cell.idents) <- levels(comb_final)
comb_final <- RenameIdents(object = comb_final, cell.idents)
comb_final$subtype <- comb_final@active.ident
comb_final <- subset(comb_final, idents = 'LQ', invert=TRUE)

#Create a New cell Class
final_obj <- comb_final
final_obj$Cell_Class <- final_obj$subtype
final_obj$Cell_Class <- as.character(final_obj$Cell_Class)
final_obj$Cell_Class[final_obj$Cell_Class %in% c('Cil', 'Basal', 'Sec', 'Gob', 'AT1', 'AT2', 'DMBT1+_AT2', "AT2_EC_DB")] <- paste0('Epithelial')
final_obj$Cell_Class[final_obj$Cell_Class %in% c('TkC', 'ThC', 'cDC', 'Mono/Macs', 'AM', 'Mast', 'GB', 'BC', 'PMN',"Perivasc_macs", "CXCL+_Immune")] <- paste0('Immune')
final_obj$Cell_Class[final_obj$Cell_Class %in% c('APOD+_FB', 'PLA2G2A+_FB','FB', 'SMC','DES+_SMC', 'Peri' )] <- paste0('Mesenchymal')
final_obj$Cell_Class[final_obj$Cell_Class %in% c('gCap','FCN3+_EC','aCap','LEC')] <- paste0('Endothelial')
saveRDS(comb_final, 'comb_final_annotated.rds')
levels(factor(final_obj$Cell_Class))
################# COPD, HC, CF Object re=clustering ############################
final_obj <- readRDS("~/Downloads/HC_COPD_CF.rds")
final_obj@active.ident <- final_obj$celltypes
final_obj <- RenameIdents(final_obj, 'CIL'='Cil') 
final_obj <- subset(final_obj, idents = c('LQ','AT2_EC_DB', 'LQ+'), invert=TRUE)
final_obj$celltypes <- final_obj@active.ident
final_obj$Cell_Class <- final_obj$celltypes
final_obj$Cell_Class <- as.character(final_obj$Cell_Class )
#Create a New cell Class
final_obj$Cell_Class[final_obj$Cell_Class %in% c('Cil', 'Bas', 'Sec', 'Gob', 'AT1', 'AT2', 'DMBT1+_AT2')] <- paste0('Epithelial')
final_obj$Cell_Class[final_obj$Cell_Class %in% c('TkC', 'ThC', 'S100B+', 'Mono/Macs', 'AM', 'Mast', 'GB', 'BC', 'PMN')] <- paste0('Immune')
final_obj$Cell_Class[final_obj$Cell_Class %in% c('STEAP4+_Peri','APOD+_FB', 'PLA2G2A+_FB','FB', 'SMC','DES+_SMC', 'Peri' )] <- paste0('Mesenchymal')
final_obj$Cell_Class[final_obj$Cell_Class %in% c('gCap','FCN+_EC','aCap','LEC')] <- paste0('Endothelial')
levels(factor(final_obj$Cell_Class))
write.csv(final_table, 'COPD_CF_HC_fractions.csv')
saveRDS(final_obj,"HC_COPD_CF_v2.rds" )

DimPlot(final_obj, label = TRUE, order = TRUE, group.by = 'Cell_Class')
FeaturePlot(final_obj, c("KIT"), cols=c("grey94","darkblue"), order=T)
VlnPlot(final_obj, c("KIT"))
##Re-cluster Cil, Giv, Sec and Bas
Ep_cluster <- subset(final_obj, subset =celltypes %in% c('Cil', 'Gob', 'Sec', 'Bas'))
Ep_cluster <- Ep_cluster %>%  
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

Ep_cluster_genes <- Ep_cluster %>% FindAllMarkers(only.pos = TRUE, recorrect_umi=FALSE)
write.csv(Ep_cluster_genes, 'ep_cluster_genes.csv')
VlnPlot(Ep_cluster, features = 'nFeature_Xenium', pt.size=0)
VlnPlot(Ep_cluster, features = 'MYH11', pt.size=0)
table(Ep_cluster$Group, Ep_cluster$seurat_clusters)
#Re-label
FeaturePlot(Ep_cluster, c('KRT5'), cols = c('grey', 'red'), order = TRUE)
cell.idents <- c('Cil', 'Cil', 'Cil', 'Gob', 'Bas', 'Bas', 'Sec_DMBT1', '?AT2', 
                 'Sec_db', 'Gob', "Bas", 'LQ_Bas', 'DB_Bas_TC', 'DB_Bas_SMC', 'DB_Bas_Mast', 
                 'DB_Gob+PMN')
Ep_cluster@active.ident <- Ep_cluster$seurat_clusters
names(cell.idents) <- levels(Ep_cluster)
Ep_cluster <- RenameIdents(object = Ep_cluster, cell.idents)
Ep_cluster$subtype <- Ep_cluster@active.ident
#Table (cells vs condition)
Ep_table <- table(Ep_cluster$subtype, Ep_cluster$Group)
Ep_table[,1] <- Ep_table[,1]/sum(Ep_table[,1])*100
Ep_table[,2] <- Ep_table[,2]/sum(Ep_table[,2])*100
Ep_table[,3] <- Ep_table[,3]/sum(Ep_table[,3])*100
Ep_table
write.csv(Ep_table, 'Ep_cell_fractions.csv')
##Re-cluster Mono/Macs, S100B+, AM, PMN
immune_cluster <- subset(final_obj, subset =celltypes %in% c('Mono/Macs', 'PMN', 'AM', 'S100B+'))
immune_cluster <- immune_cluster %>%  
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)
DimPlot(immune_cluster, label=TRUE)
immune_cluster_genes <- immune_cluster %>% FindAllMarkers(only.pos = TRUE, recorrect_umi=FALSE)
write.csv(immune_cluster_genes, 'immune_cluster_genes.csv')
VlnPlot(immune_cluster, features = 'CHIT1', pt.size=0)
#Re-label
cell.idents <- c('Macrophages', 'Perivascular Mac', 'Mono', 'AM', 'PI3+ PMN', 'AM',
                 'PMN', 'cDC', 'EC_db', 'FB_db', 'DC2', 'AT2_db', 'PMN', 'LQ_perivasc mac',
                 'LQ_Pi3_hi PMN', 'GB_db', 'PMN', 'SMC_db', 'PMN')
immune_cluster@active.ident <- immune_cluster$seurat_clusters
names(cell.idents) <- levels(immune_cluster)
immune_cluster <- RenameIdents(object = immune_cluster, cell.idents)
immune_cluster$subtype <- immune_cluster@active.ident
#Table (cells vs condition)
immune_table <- table(immune_cluster$subtype, immune_cluster$Group)
immune_table[,1] <- immune_table[,1]/sum(immune_table[,1])*100
immune_table[,2] <- immune_table[,2]/sum(immune_table[,2])*100
immune_table[,3] <- immune_table[,3]/sum(immune_table[,3])*100
immune_table
write.csv(immune_table, 'immune_cell_fractions.csv')
FeaturePlot(immune_cluster,c('MMP12'),  cols = c('grey94', 'red'), order=T)


saveRDS(immune_cluster, 'immune_cluster.rds')
saveRDS(Ep_cluster, 'ep_cluster.rds')

###############BarPlots#########################################################
library(RColorBrewer)
coul <- brewer.pal(5, "Set3")
setwd("~/Spatial_Data/BarPlots")
pdf('COPD_HC_CF_Counts_per_Group.pdf', 8,8)
barplot(final_table,xlim=c(0, ncol(final_table) + 2),col=brewer.pal(nrow(final_table), "Paired"), ylab = 'Fraction of Cells', main='Celltypes Across Lung Samples', xlab='Experimental Group',legend.text=TRUE,
        args.legend=list( x=ncol(final_table) + 2,y=max(colSums(final_table)),bty = "n"))
dev.off()

setwd("/media/haider/Files/Spatial_Data/BarPlots_2")
final_obj <- readRDS("/media/haider/Files/Spatial_Data/annotated_comb_final.rds")
final_obj$Group <- factor(x = final_obj$Group, levels = c("HC", 'COPD', 'CF', 'Asthma', 'BOS'))
#Immune Cells
immune_cluster <- subset(final_obj, subset = Cell_Class %in% c('Immune'))
immune_table <- table(immune_cluster@active.ident, immune_cluster$Group)
immune_table[,1] <- immune_table[,1]/sum(immune_table[,1])*100
immune_table[,2] <- immune_table[,2]/sum(immune_table[,2])*100
immune_table[,3] <- immune_table[,3]/sum(immune_table[,3])*100
immune_table[,4] <- immune_table[,4]/sum(immune_table[,4])*100
immune_table[,5] <- immune_table[,5]/sum(immune_table[,5])*100
immune_table
#Epithelial Cells
ep_cluster <- subset(final_obj, subset = Cell_Class %in% c("Epithelial"))
ep_table <- table(ep_cluster@active.ident, ep_cluster$Group)
ep_table[,1] <- ep_table[,1]/sum(ep_table[,1])*100
ep_table[,2] <- ep_table[,2]/sum(ep_table[,2])*100
ep_table[,3] <- ep_table[,3]/sum(ep_table[,3])*100
ep_table[,4] <- ep_table[,4]/sum(ep_table[,4])*100
ep_table[,5] <- ep_table[,5]/sum(ep_table[,5])*100
ep_table
#Endothelial Cells
end_cluster <- subset(final_obj, subset = Cell_Class %in% c("Endothelial"))
end_table <- table(end_cluster@active.ident, end_cluster$Group)
end_table[,1] <- end_table[,1]/sum(end_table[,1])*100
end_table[,2] <- end_table[,2]/sum(end_table[,2])*100
end_table[,3] <- end_table[,3]/sum(end_table[,3])*100
end_table[,4] <- end_table[,4]/sum(end_table[,4])*100
end_table[,5] <- end_table[,5]/sum(end_table[,5])*100
end_table
#Mesenchymal Cells
mes_cluster <- subset(final_obj, subset = Cell_Class %in% c("Mesenchymal"))
mes_table <- table(mes_cluster@active.ident, mes_cluster$Group)
mes_table[,1] <- mes_table[,1]/sum(mes_table[,1])*100
mes_table[,2] <- mes_table[,2]/sum(mes_table[,2])*100
mes_table[,3] <- mes_table[,3]/sum(mes_table[,3])*100
mes_table[,4] <- mes_table[,4]/sum(mes_table[,4])*100
mes_table[,5] <- mes_table[,5]/sum(mes_table[,5])*100
mes_table

pdf("COPD_HC_CF_BOS_Asthma_SubCluster_Counts_per_Group.pdf", 15, 15)
barplot(immune_table,xlim=c(0, ncol(immune_table) + 2),col=brewer.pal(nrow(immune_table), "Paired"), ylab = 'Fraction of Cells', main='Immune Cells', xlab='Experimental Group',legend.text=TRUE,
        args.legend=list( x=ncol(immune_table) + 2,y=max(colSums(immune_table)),bty = "n"))

barplot(end_table,xlim=c(0, ncol(end_table) + 2),col=brewer.pal(nrow(end_table), "Paired"), ylab = 'Fraction of Cells', main='Endothelial Cells', xlab='Experimental Group',legend.text=TRUE,
        args.legend=list( x=ncol(end_table) + 2,y=max(colSums(end_table)),bty = "n"))

barplot(ep_table,xlim=c(0, ncol(ep_table) + 2),col=brewer.pal(nrow(ep_table), "Paired"), ylab = 'Fraction of Cells', main='Epithelial Cells', xlab='Experimental Group',legend.text=TRUE,
        args.legend=list( x=ncol(ep_table) + 2,y=max(colSums(ep_table)),bty = "n"))

barplot(mes_table,xlim=c(0, ncol(mes_table) + 2),col=brewer.pal(nrow(mes_table), "Paired"), ylab = 'Fraction of Cells', main='Mesenchymal Cells', xlab='Experimental Group',legend.text=TRUE,
        args.legend=list( x=ncol(mes_table) + 2,y=max(colSums(mes_table)),bty = "n"))
dev.off()
#saveRDS(final_obj, 'annotated_comb_final.rds')
write.csv(immune_table, 'immune_cell_fractions_final.csv')
write.csv(end_table, 'endothelial_cell_fractions_final.csv')
write.csv(ep_table, 'epithelial_cell_fractions_final.csv')
write.csv(mes_table, 'mesenchymal_cell_fractions_final.csv')



#################Create Boxplots for Disease Comparison########################
library(ggplot2)
library(ggrepel)
final_table2 <- table(final_obj@active.ident, final_obj$ID)
final_table2[,1] <- final_table2[,1]/sum(final_table2[,1])*100
final_table2[,2] <- final_table2[,2]/sum(final_table2[,2])*100
final_table2[,3] <- final_table2[,3]/sum(final_table2[,3])*100
final_table2[,4] <- final_table2[,4]/sum(final_table2[,4])*100
final_table2[,5] <- final_table2[,5]/sum(final_table2[,5])*100
final_table2[,6] <- final_table2[,6]/sum(final_table2[,6])*100
final_table2[,7] <- final_table2[,7]/sum(final_table2[,7])*100
final_table2[,8] <- final_table2[,8]/sum(final_table2[,8])*100
final_table2[,9] <- final_table2[,9]/sum(final_table2[,9])*100
final_table2[,10] <- final_table2[,10]/sum(final_table2[,10])*100
final_table2[,11] <- final_table2[,11]/sum(final_table2[,11])*100
final_table2[,12] <- final_table2[,12]/sum(final_table2[,12])*100
final_table2[,13] <- final_table2[,13]/sum(final_table2[,13])*100
final_table2[,14] <- final_table2[,14]/sum(final_table2[,14])*100
final_table2[,15] <- final_table2[,15]/sum(final_table2[,15])*100
final_table2[,16] <- final_table2[,16]/sum(final_table2[,16])*100
final_table2[,17] <- final_table2[,17]/sum(final_table2[,17])*100
final_table2[,18] <- final_table2[,18]/sum(final_table2[,18])*100
final_table2[,19] <- final_table2[,19]/sum(final_table2[,19])*100
final_table2[,20] <- final_table2[,20]/sum(final_table2[,20])*100
final_table2[,21] <- final_table2[,21]/sum(final_table2[,21])*100
final_table2[,22] <- final_table2[,22]/sum(final_table2[,22])*100
final_table2 <- as.data.frame(final_table2)

final_table2$Var2 <- as.character(final_table2$Var2)
final_table2[final_table2$Var2 %in% c('BOS11A', 'BOS12D', 'BOS13E'),]$Var2 <- 'BOS'
final_table2[final_table2$Var2 %in% c('MTS21A', 'MTS48B', 'MTS27E'),]$Var2 <- 'Asthma'
final_table2[final_table2$Var2 %in% c('ID6','MTS7D', 'MTS12D', 'MTS29E', 'MTS38D', 'MTS43D'),]$Var2 <- 'HC'
final_table2[final_table2$Var2 %in% c('C48A','C50C','C51E' ,'C52E', 'C56E'),]$Var2 <- 'COPD'
final_table2[final_table2$Var2 %in% c('CF19', 'CF20', 'CF27', 'CF35', 'CF40'),]$Var2 <- 'CF'
final_table2$Var2 <- as.factor(final_table2$Var2)
final_table2$Var2 <- factor(final_table2$Var2, levels = c('HC',  "Asthma" ,"BOS"   , "CF" ,    "COPD" ))

immune_list <- c('TkC', 'ThC', 'cDC', 'AM', 'Mast', 'GB', 'BC', 'PMN',"Perivasc_macs", "CXCL+_Immune", 'Mono/Macs')
alveolar_list <- c('AT1', 'AT2', 'AT2_EC_DB', 'DMBT1+_AT2', 'aCap', 'gCap','FCN3+_EC')
vascular_list <- c('Basal', 'Cil', 'Sec', 'Gob', 'SMC', 'DES+_SMC', 'PLA2G2A+_FB')
mes_list <- c('APOD+_FB','FB','Peri', 'LEC')

immune_freq_table <- final_table2[final_table2$Var1 %in% immune_list,]
vascular_freq_table <- final_table2[final_table2$Var1 %in% vascular_list,]
write.csv(immune_freq_table, "~/Squidpy_Spatial_analysis/Boxplots/Frequency_table/immune.csv")
write.csv(vascular_freq_table, "~/Squidpy_Spatial_analysis/Boxplots/Frequency_table/vascular.csv")

whisker <- function(sample_list, obj){
  for (i in 1:length(sample_list)){
    print(ggplot(obj[obj$Var1 %in% c(sample_list[i]),], aes(x=`Var2`, y=Freq)) +
            geom_boxplot(outlier.shape = NA)+ geom_jitter(width = 0.2)+ggtitle(sample_list[i]))
}}

setwd("~/Squidpy_Spatial_analysis/Boxplots/immune")
pdf('immune_list_box_allcond.pdf',5,5)
whisker(immune_list, final_table2)
dev.off()s
whisker(alveolar_list, final_table2)
whisker(vascular_list, final_table2)
whisker(mes_list, final_table2)

#Immune Cells
#immune_cluster <- subset(final_obj, subset = Cell_Class %in% c('Immune'))
immune_table <- table(immune_cluster@active.ident, immune_cluster$ID)

cell_fractions <- function(table_obj){
  for (i in 1:22){
    table_obj[,i] <- table_obj[,i]/sum(table_obj[,i])*100}
  table_obj <- as.data.frame(table_obj)
  table_obj$Var2 <- as.character(table_obj$Var2)
  table_obj[table_obj$Var2 %in% c('BOS11A', 'BOS12D', 'BOS13E'),]$Var2 <- 'BOS'
  table_obj[table_obj$Var2 %in% c('MTS21A', 'MTS48B', 'MTS27E'),]$Var2 <- 'Asthma'
  table_obj[table_obj$Var2 %in% c('ID6','MTS7D', 'MTS12D', 'MTS29E', 'MTS38D', 'MTS43D'),]$Var2 <- 'HC'
  table_obj[table_obj$Var2 %in% c('C48A','C50C','C51E' ,'C52E', 'C56E'),]$Var2 <- 'COPD'
  table_obj[table_obj$Var2 %in% c('CF19', 'CF20', 'CF27', 'CF35', 'CF40'),]$Var2 <- 'CF'
  table_obj$Var2 <- as.factor(table_obj$Var2)
  return(table_obj)}


# Normalize based on total basal cells
final_table2 <- table(final_obj@active.ident, final_obj$ID)

cell_fractions_bas <- function(table_obj){
  for (i in 1:22){
    #bas_percent <- table_obj[10,i]/sum(table_obj[,i])*100
    cil_percent <- table_obj[5,i]/sum(table_obj[,i])*100
    table_obj[,i] <- (table_obj[,i]/sum(table_obj[,i])*100)/cil_percent}
  table_obj <- as.data.frame(table_obj)
  table_obj$Var2 <- as.character(table_obj$Var2)
  table_obj[table_obj$Var2 %in% c('BOS11A', 'BOS12D', 'BOS13E'),]$Var2 <- 'BOS'
  table_obj[table_obj$Var2 %in% c('MTS21A', 'MTS48B', 'MTS27E'),]$Var2 <- 'Asthma'
  table_obj[table_obj$Var2 %in% c('ID6','MTS7D', 'MTS12D', 'MTS29E', 'MTS38D', 'MTS43D'),]$Var2 <- 'HC'
  table_obj[table_obj$Var2 %in% c('C48A','C50C','C51E' ,'C52E', 'C56E'),]$Var2 <- 'COPD'
  table_obj[table_obj$Var2 %in% c('CF19', 'CF20', 'CF27', 'CF35', 'CF40'),]$Var2 <- 'CF'
  table_obj$Var2 <- as.factor(table_obj$Var2)
  table_obj$Var2 <- factor(table_obj$Var2, levels = c('HC', 'Asthma', 'BOS', 'CF', 'COPD'))
  return(table_obj)}
final_table2 <- cell_fractions_bas(final_table2)
whisker(immune_list, final_table2)
whisker(vascular_list, final_table2)

################################################################################
levels(factor(final_obj@meta.data$ID))

##Tables of cell types by sample ID (immune)
immune_table2 <- table(immune_cluster@active.ident, immune_cluster$ID)
immune_table2[,1] <- immune_table2[,1]/sum(immune_table2[,1])*100
immune_table2[,2] <- immune_table2[,2]/sum(immune_table2[,2])*100
immune_table2[,3] <- immune_table2[,3]/sum(immune_table2[,3])*100
immune_table2[,4] <- immune_table2[,4]/sum(immune_table2[,4])*100
immune_table2[,5] <- immune_table2[,5]/sum(immune_table2[,5])*100
immune_table2[,6] <- immune_table2[,6]/sum(immune_table2[,6])*100
immune_table2[,7] <- immune_table2[,7]/sum(immune_table2[,7])*100
immune_table2[,8] <- immune_table2[,8]/sum(immune_table2[,8])*100
immune_table2[,9] <- immune_table2[,9]/sum(immune_table2[,9])*100
immune_table2[,10] <- immune_table2[,10]/sum(immune_table2[,10])*100
immune_table2[,11] <- immune_table2[,11]/sum(immune_table2[,11])*100
immune_table2[,12] <- immune_table2[,12]/sum(immune_table2[,12])*100
immune_table2[,13] <- immune_table2[,13]/sum(immune_table2[,13])*100
immune_table2[,14] <- immune_table2[,14]/sum(immune_table2[,14])*100
immune_table2[,15] <- immune_table2[,15]/sum(immune_table2[,15])*100
immune_table2[,16] <- immune_table2[,16]/sum(immune_table2[,16])*100
immune_table2

#Endothelial
end_table2 <- table(end_cluster@active.ident, end_cluster$ID)
end_table2[,1] <- end_table2[,1]/sum(end_table2[,1])*100
end_table2[,2] <- end_table2[,2]/sum(end_table2[,2])*100
end_table2[,3] <- end_table2[,3]/sum(end_table2[,3])*100
end_table2[,4] <- end_table2[,4]/sum(end_table2[,4])*100
end_table2[,5] <- end_table2[,5]/sum(end_table2[,5])*100
end_table2[,6] <- end_table2[,6]/sum(end_table2[,6])*100
end_table2[,7] <- end_table2[,7]/sum(end_table2[,7])*100
end_table2[,8] <- end_table2[,8]/sum(end_table2[,8])*100
end_table2[,9] <- end_table2[,9]/sum(end_table2[,9])*100
end_table2[,10] <- end_table2[,10]/sum(end_table2[,10])*100
end_table2[,11] <- end_table2[,11]/sum(end_table2[,11])*100
end_table2[,12] <- end_table2[,12]/sum(end_table2[,12])*100
end_table2[,13] <- end_table2[,13]/sum(end_table2[,13])*100
end_table2[,14] <- end_table2[,14]/sum(end_table2[,14])*100
end_table2[,15] <- end_table2[,15]/sum(end_table2[,15])*100
end_table2[,16] <- end_table2[,16]/sum(end_table2[,16])*100
end_table2

#Epithelial
ep_table2 <- table(ep_cluster@active.ident, ep_cluster$ID)
ep_table2[,1] <- ep_table2[,1]/sum(ep_table2[,1])*100
ep_table2[,2] <- ep_table2[,2]/sum(ep_table2[,2])*100
ep_table2[,3] <- ep_table2[,3]/sum(ep_table2[,3])*100
ep_table2[,4] <- ep_table2[,4]/sum(ep_table2[,4])*100
ep_table2[,5] <- ep_table2[,5]/sum(ep_table2[,5])*100
ep_table2[,6] <- ep_table2[,6]/sum(ep_table2[,6])*100
ep_table2[,7] <- ep_table2[,7]/sum(ep_table2[,7])*100
ep_table2[,8] <- ep_table2[,8]/sum(ep_table2[,8])*100
ep_table2[,9] <- ep_table2[,9]/sum(ep_table2[,9])*100
ep_table2[,10] <- ep_table2[,10]/sum(ep_table2[,10])*100
ep_table2[,11] <- ep_table2[,11]/sum(ep_table2[,11])*100
ep_table2[,12] <- ep_table2[,12]/sum(ep_table2[,12])*100
ep_table2[,13] <- ep_table2[,13]/sum(ep_table2[,13])*100
ep_table2[,14] <- ep_table2[,14]/sum(ep_table2[,14])*100
ep_table2[,15] <- ep_table2[,15]/sum(ep_table2[,15])*100
ep_table2[,16] <- ep_table2[,16]/sum(ep_table2[,16])*100
ep_table2

#Mesenchymal
mes_table2 <- table(mes_cluster@active.ident, mes_cluster$ID)
mes_table2[,1] <- mes_table2[,1]/sum(mes_table2[,1])*100
mes_table2[,2] <- mes_table2[,2]/sum(mes_table2[,2])*100
mes_table2[,3] <- mes_table2[,3]/sum(mes_table2[,3])*100
mes_table2[,4] <- mes_table2[,4]/sum(mes_table2[,4])*100
mes_table2[,5] <- mes_table2[,5]/sum(mes_table2[,5])*100
mes_table2[,6] <- mes_table2[,6]/sum(mes_table2[,6])*100
mes_table2[,7] <- mes_table2[,7]/sum(mes_table2[,7])*100
mes_table2[,8] <- mes_table2[,8]/sum(mes_table2[,8])*100
mes_table2[,9] <- mes_table2[,9]/sum(mes_table2[,9])*100
mes_table2[,10] <- mes_table2[,10]/sum(mes_table2[,10])*100
mes_table2[,11] <- mes_table2[,11]/sum(mes_table2[,11])*100
mes_table2[,12] <- mes_table2[,12]/sum(mes_table2[,12])*100
mes_table2[,13] <- mes_table2[,13]/sum(mes_table2[,13])*100
mes_table2[,14] <- mes_table2[,14]/sum(mes_table2[,14])*100
mes_table2[,15] <- mes_table2[,15]/sum(mes_table2[,15])*100
mes_table2[,16] <- mes_table2[,16]/sum(mes_table2[,16])*100
mes_table2

write.csv(immune_table2, 'immune_ID_fractions.csv')
write.csv(end_table2, 'end_ID_fractions.csv')
write.csv(ep_table2, 'ep_ID_fractions.csv')
write.csv(mes_table2, 'mes_ID_fractions.csv')

final_obj <- readRDS('~/Spatial_Data/annotated_comb_final.rds')
seurat_cord <-final_obj@meta.data[c('Cell_Class', 'subtype', 'seurat_clusters', 'ID', 'Group')]
seurat_cord$X <- sapply(strsplit(as.character(rownames(seurat_cord)), "_"), `[`, 2)
seurat_cord <- seurat_cord[!duplicated(seurat_cord$X),]
write.csv(seurat_cord, "~/Squidpy_Spatial_analysis/seurat_cord.csv")

# Generate DEG list based on Disease Group
final_obj@active.ident <- factor(final_obj$Group)
setwd("~/Spatial_Data/immune_ep_DEGS")
Cell_list <- c('PMN','AM', 'Mast', 'cDC', 'GB', 'BC', 'Cil' ,'Basal', 'Sec')
for (i in 1:length(Cell_list)){
  obj <- subset(final_obj, subset = subtype %in% Cell_list[i])s
  file <- obj %>% FindAllMarkers(only.pos = TRUE,  recorrect_umi=FALSE)
  write.csv(file, paste0(Cell_list[i], '_DG.csv'))
}
# Mono/Macs
Mono <- subset(final_obj, subset = subtype %in% c('Mono/Macs'))
Mono_genes <- Mono %>% FindAllMarkers(only.pos = TRUE, recorrect_umi=FALSE)
write.csv(Mono_genes, 'Mono_DG.csv')
# T cells
T_cell <- subset(final_obj, subset = subtype %in% c('TkC', 'ThC'))
T_genes <- T_cell %>% FindAllMarkers(only.pos = TRUE, recorrect_umi=FALSE)
write.csv(T_genes, 'Tcell_DG.csv')


