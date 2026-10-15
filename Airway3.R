library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
library(anndata)
library(ggrepel)
library(scCustomize)
library(anndata)
plan(multisession(workers = 20))
options(future.globals.maxSize = 8000 * 1024^3)
`%ni%` <- Negate(`%in%`)
library(reticulate)
library(SeuratExtend)
library(SeuratDisk)
library(ggrepel)
library(circlize)

library(reticulate)
path_to_python <- "/home/haider/miniconda3/bin/python"
use_python(path_to_python)
################################################################################p
# Load complete spatial object
final_obj <- readRDS("/media/haider/Files/Spatial_Data/annotated_comb_final.rds")
final_obj <- UpdateSeuratObject(final_obj)
# Extract the raw expression matrix from an existing Seurat object
expression_matrix <- GetAssayData(object = ep_obj, assay = "SCT", slot = "counts")

# Create a new Seurat object from the extracted matrix
new_ep_obj <- CreateSeuratObject(counts = expression_matrix, project = "NewProject")
new_ep_obj[['SCT']] <- ep_obj[['SCT']]
new_ep_obj @meta.data <- ep_obj @meta.data
new_ep_obj @reductions$umap <- ep_obj @reductions$umap
# Remove low quality clusters 
new_ep_obj  <- subset(new_ep_obj , subset=subtype %ni% c('LQ',  "AT2_EC_DB" ))
new_ep_obj  <- subset(new_ep_obj , subset= ID %ni% c('MTS43D', 'MTS38D'))


ep_obj <- readRDS("/media/haider/Files/Squidpy_Spatial_analysis/ep_obj.rds")
#
Mast_all_layers <- subset(new_final_obj, subset=subtype %in% c('Mast' ))
acap_all_layers <- subset(new_final_obj, subset=subtype %in% c('aCap'))

# Remove low quality clusters 
final_obj <- subset(final_obj, subset=subtype %ni% c('LQ',  "AT2_EC_DB" ))
final_obj <- subset(final_obj, subset= ID %ni% c('MTS43D', 'MTS38D'))

# Add layers
data <- anndata::read_h5ad("/media/haider/Files/Squidpy_Spatial_analysis/Airway_subclass_obj")
ep <- data$obs[data$obs$Layers %in% 'Airway_Ep',]
sub <- data$obs[data$obs$Layers %in% 'Sub_Ep',]
smg <- data$obs[data$obs$Layers %in% 'SMG',]

final_obj$Layers <- sapply(strsplit(as.character(colnames(final_obj)), "_"), `[`, 2)
final_obj$Layers[final_obj$Layers %in% ep$cell_ids] <- 'Airway_Ep'
final_obj$Layers[final_obj$Layers %in% sub$cell_ids] <- 'Sub_Ep'
final_obj$Layers[final_obj$Layers %in% smg$cell_ids] <- 'SMG'
final_obj$Layers[final_obj$Layers %ni% c('Airway_Ep', 'Sub_Ep', 'SMG')] <- 'Non-Airway'

# Epithelial and Sub-Epithelial Layer
final_obj <- readRDS("/media/haider/Files/Squidpy_Spatial_analysis/final_annot_layers.rds")
ep_sub <- subset(final_obj, subset=Layers %in% c('Airway_Ep', 'Sub_Ep'))
sub <- subset(final_obj, subset=Layers %in% c('Sub_Ep'))


# Name NK cells
#subcluster the TkC cells
TkC <- FindSubCluster(final_obj, cluster = 'TkC', graph.name = 'SCT_nn', resolution = 0.3)
NK_cell_ids <- rownames(TkC@meta.data[TkC$sub.cluster %in% 'TkC_4',])
final_obj$subtype <- as.character(final_obj$subtype)
final_obj$subtype[rownames(final_obj@meta.data) %in% NK_cell_ids] <- paste0('NK')
final_obj$subtype <- as.factor(final_obj$subtype)
final_obj@active.ident <- final_obj$subtype
final_obj<- RenameIdents(final_obj , 'TkC'='CD8 Tcells')
final_obj<- RenameIdents(final_obj , 'ThC'='CD4 Tcells')

# Featrureplot across disease types for cxcl immune cells
pdf('/media/haider/Files/Squidpy_Spatial_analysis/CXCL_IMMUNE_BY_DISEASE/Featureplot_across_disease_allcxcl.pdf',8,8)
FeaturePlot(final_obj, features = c('CXCL9'), cols = c('grey94', 'darkblue'), split.by = 'Group')
dev.off()


# GATA3 positive cells subset
gata3_sub_ep <- subset(ep_sub, subset=subtype %in% c('ThC', 'TkC', 'CXCL+_Immune'))

# CXCL Immune cluster subset 
cxcl_all <- subset(final_obj, subset=subtype %in% c('CXCL+_Immune'))
Idents(cxcl_all) <- cxcl_all$Group
DoHeatmap(AverageExpression(cxcl_all, return.seurat = T, assays ='SCT' ), features = c('CXCL9', 'CXCL10', 'CXCL13', 'CD3E', 'EGFR', 'PDGFRA', 'AGER', 'ACTA', 
                                 'CD34', 'CD19', 'FOXP3','GATA3','KIT'), draw.lines = F)+ scale_fill_gradientn(colors=c('blue', 'white', 'red'))

# Subset ILC2 cells
gata3 <- subset(final_obj, subset=subtype %in% c('ThC', 'TkC', 'CXCL+_Immune'))
gata3_sub <- gata3 %>%  FindVariableFeatures() %>%
  ScaleData() %>% 
  RunPCA(verbose = FALSE) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)
Immune_sub_markers <- FindAllMarkers(Immune_sub, only.pos = T, recorrect_umi=FALSE)



# Ratio of Mast cells (Ep vs Sub Ep across Disease )
mast_ep_sub <- subset(ep_sub, subset=subtype %in% 'Mast')
mast_ep_sub_df <- table(mast_ep_sub$Layers, mast_ep_sub$Group)
cum_vals <- apply(mast_ep_sub_df, 2, cumsum)

# 🔹 Give space on right + allow drawing outside plot
par(mar = c(5, 4, 4, 8), xpd = TRUE)

bp <- barplot(mast_ep_sub_df, 
              main = "Mast Cell Ratio (Epithelial vs Sub-Epithelial) by Disease Group",
              xlab = "Disease Group", 
              ylab = "Mast Cell Frequency", 
              border = NA,
              col = c('blue', 'darkgreen')
)

legend("topright",
       inset = c(-0.15, 0),   # move outside
       legend = c("Airway Epithelium", "Sub-Epithelial"),
       fill = c("blue", "darkgreen"),
       title = "Layer",
       bty = "n",
       cex = 0.9
)

text(
  x = rep(bp, each = nrow(mast_ep_sub_df)),
  y = as.vector(cum_vals - mast_ep_sub_df / 2),
  labels = as.vector(mast_ep_sub_df),
  cex = 1,
  font = 2,
  col = 'white'
)


# Convert to proportions by column (within each disease group)
mast_prop <- prop.table(mast_ep_sub_df, margin = 2)

# Add margin space on right for legend
par(mar = c(5, 4, 4, 8), xpd = TRUE)

bp <- barplot(mast_prop,
              main = "Mast Cell Ratio (Epithelial vs Sub-Epithelial) by Disease Group",
              xlab = "Disease Group",
              ylab = "Proportion of Mast Cells",
              col = c("blue", "darkgreen"),
              ylim = c(0, 1),
              border = NA
)

legend("topright",
       inset = c(-0.15, 0),
       legend = c("Airway Epithelium", "Sub-Epithelial"),
       fill = c("blue", "darkgreen"),
       title = "Layer",
       bty = "n"
)

# Add proportion labels inside bars
cum_vals <- apply(mast_prop, 2, cumsum)

text(
  x = rep(bp, each = nrow(mast_prop)),
  y = as.vector(cum_vals - mast_prop / 2),
  labels = round(as.vector(mast_prop), 2),
  col = "white",
  font = 2
)


# TMEM16 Splice variants highlighted on UMAP
set.seed(1)
cluster_rand <- rand_color(length(levels(factor(final_obj$subtype))))
cf_hc_subset <- subset(final_obj, Group %in% c('HC', 'CF'))

pdf('CF_HC_DimPlot.pdf', 8,8)
DimPlot(cf_hc_subset, label=T, order=T, raster = F, cols = cluster_rand)
dev.off()
pdf('CF_HC_ANO1_Featureplot.pdf',8,8)
FeaturePlot(cf_hc_subset,features = c('ANO1'), split.by = 'Group', raster = F)
dev.off()

saveRDS(final_obj, 'final_annot_layers.rds')

# HC vs CF (CXCL+_Immune and BC)
CF_CXCL_ep_sub <- subset(ep_sub, subset=Group %in% c('CF', 'HC')) %>%  subset(subtype %in% c("CXCL+_Immune"))
CF_BC_ep_sub <- subset(ep_sub, subset=Group %in% c('CF', 'HC')) %>%  subset(subtype %in% c("BC"))

Bos_ep_sub <- subset(ep_sub, subset=Group %in% c('BOS', 'HC')) %>% subset( subset=subtype %in% c("Mono/Macs"))


Idents(CF_CXCL_ep_sub)<-CF_ep_sub$Group
CF_CXCL_ep_sub.df <- FindMarkers(CF_CXCL_ep_sub, ident.1 = 'CF', ident.2 = 'HC',recorrect_umi=FALSE)
CF_CXCL_ep_sub.df$gene <- rownames(CF_CXCL_ep_sub.df)

# Immune subset (cDC, Mono/macs, AM, PMN, Perivasc_Macs)
Immune_sub <- subset(final_obj, subset=subtype %in% c('cDC', 'Perivasc_macs', 'PMN', 'AM', 'Mono/Macs'))
Immune_sub_int <- Immune_sub %>%
  RunPCA(npcs = 30, verbose = FALSE)
ElbowPlot(Immune_sub_int)
# one-liner to run Integration
Immune_sub_int = IntegrateLayers(object =Immune_sub_int, method = "HarmonyIntegration", normalization.method = "SCT", verbose = T, new.reduction = "harmony", orig.reduction = "pca")
Immune_sub_int <- Immune_sub_int %>%  
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.4, verbose = FALSE)

Immune_sub_int_markers <- FindAllMarkers(Immune_sub_int, only.pos=T, max.cells.per.ident=500, recorrect_umi=FALSE)
write.csv(Immune_sub_int_markers, file="Immune_sub_int_markers.csv")

Immune_sub <- subset(final_obj, subset=subtype %in% c('cDC', 'Perivasc_macs', 'PMN', 'AM', 'Mono/Macs'))
Immune_sub <- Immune_sub %>%  
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.4, verbose = FALSE)
Immune_sub_markers <- FindAllMarkers(Immune_sub, only.pos = T, recorrect_umi=FALSE)
write.csv(Immune_sub_markers, "~/Squidpy_Spatial_analysis/Mono_mac_recluster/Immune_mono_markers.csv")

Idents(Immune_sub_int) <- "seurat_clusters"
DimPlot(Immune_sub_int, label=TRUE)
saveRDS(Immune_sub, 'macs_mono_cluster.rds')
FeaturePlot(Immune_sub_int, c("SLC1A3"), cols=c("grey94","darkblue"), order=T)
DimPlot(final_obj, label=T) + NoLegend()

# Remove Non-Specific clusters from the immune object
Immune_sub_int_spec <- subset(Immune_sub_int, subset=seurat_clusters %ni% c(22,23,25,26,27,28,30,32,33))
Immune_sub_int_spec <- RenameIdents(object = Immune_sub_int_spec, '31' = '1')
Immune_sub_int_spec <- RenameIdents(object = Immune_sub_int_spec, '24' = '3')
Immune_sub_int_spec <- RenameIdents(object = Immune_sub_int_spec, '10' = '4')
DimPlot(Immune_sub_int_spec, label=TRUE)
Immune_sub_int_spec$seurat_clusters2 <- Immune_sub_int_spec@active.ident 


FeaturePlot(final_obj, c("MPEG1"), cols=c("grey94","darkblue"), order=T)

# BoxPlots
create_table <- function(final_table2){
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
  final_table2 <- as.data.frame(final_table2)
  final_table2$Var2 <- as.character(final_table2$Var2)
  final_table2[final_table2$Var2 %in% c('BOS11A', 'BOS12D', 'BOS13E'),]$Var2 <- 'BOS'
  final_table2[final_table2$Var2 %in% c('MTS21A', 'MTS48B', 'MTS27E'),]$Var2 <- 'Asthma'
  final_table2[final_table2$Var2 %in% c('ID6','MTS7D', 'MTS12D', 'MTS29E'),]$Var2 <- 'HC'
  final_table2[final_table2$Var2 %in% c('C48A','C50C','C51E' ,'C52E', 'C56E'),]$Var2 <- 'COPD'
  final_table2[final_table2$Var2 %in% c('CF19', 'CF20', 'CF27', 'CF35', 'CF40'),]$Var2 <- 'CF'
  final_table2$Var2 <- as.factor(final_table2$Var2)
  final_table2$Var2 <- factor(final_table2$Var2, levels = c('HC', 'Asthma', 'BOS', 'CF', 'COPD'))
  return(final_table2)
}

# Immune integrated table 
Immune_sub_int@active.ident <- factor(Immune_sub_int$Group)
Immune_0_5 <- subset(Immune_sub_int, subset = seurat_clusters %in% c(0,1,2,3,4,5))
Immune_int_table <- create_table(table(Immune_0_5$seurat_clusters, Immune_0_5$ID))

# Immune integrated table with non specific gene removed
Immune_int_spec_table <- create_table(table(Immune_sub_int_spec$seurat_clusters2, Immune_sub_int_spec$ID))

Immune_sub_int_spec_ep <- subset(Immune_sub_int_spec, subset=Layers %in% c('Airway_Ep', 'Sub_Ep'))
Immune_int_spec_ep_table <- create_table(table(Immune_sub_int_spec_ep$seurat_clusters2, Immune_sub_int_spec_ep$ID))


whisker <- function(sample_list, obj){
  for (i in 1:length(sample_list)){
    print(ggplot(obj[obj$Var1 %in% c(sample_list[i]),], aes(x=`Var2`, y=Freq)) +
            geom_boxplot(outlier.shape = NA)+ geom_jitter(width = 0.2)+ggtitle(sample_list[i]))
  }}

saveRDS(Immune_sub_int_spec, 'Immune_sub_int_spec.rds')
pdf('Immune_Box_integrated_plots', 5,5)
whisker(c("0", "1",  "2" , "3", "4" , "5" , "6" , "7" , "8" , "9" , "10" ,"11", "12" ,"13", "14", "15", "16", "17", "18", "19", "20", "21", "29", "34"), Immune_int_spec_table)
dev.off()

saveRDS(Immune_sub_int_spec_ep, 'Immune_sub_int_ep_spec.rds')
pdf('Immune_Box_integrated_ep_sub_plots', 5,5)
whisker(c("0", "1",  "2" , "3", "4" , "5" , "6" , "7" , "8" , "9" , "10" ,"11", "12" ,"13", "14", "15", "16", "17", "18", "19", "20", "21", "29", "34"), Immune_int_spec_ep_table)
dev.off()


DimPlot(Immune_sub_int_spec_ep, split.by="Group", label=T)
FeaturePlot(Immune_sub_int_spec_ep, c("MS4A4A"), cols=c("grey94","darkblue"))

# Mast cells Sub-Ep vs Airway Epithelial
Idents(mast_ep_sub)<-mast_ep_sub$Layers
mast_ep_sub.df <- FindMarkers(mast_ep_sub, ident.1 =  "Airway_Ep" , ident.2 = "Sub_Ep" ,recorrect_umi=FALSE)
mast_ep_sub.df$gene <- rownames(mast_ep_sub.df)

# GATA3 cells Sub-Ep vs Airway Epithelial
Idents(gata3_sub_ep)<-gata3_sub_ep$Layers
gata3_sub_ep.df <- FindMarkers(gata3_sub_ep, ident.1 =  "Airway_Ep" , ident.2 = "Sub_Ep" ,recorrect_umi=FALSE)
gata3_sub_ep.df$gene <- rownames(gata3_sub_ep.df)


# Mast all Layers DE
Idents(Mast_all_layers )<-Mast_all_layers$Group
Mast_all_layers.df <- FindMarkers(Mast_all_layers, ident.1 = 'COPD', ident.2 = 'HC',recorrect_umi=FALSE)
Mast_all_layers.df$gene <- rownames(Mast_all_layers.df)


#acap all Layers DE
Idents(acap_all_layers)<-acap_all_layers$Group
acap_all_layers.df <- FindMarkers(acap_all_layers, ident.1 = 'COPD', ident.2 = 'HC',recorrect_umi=FALSE)
acap_all_layers.df$gene <- rownames(acap_all_layers.df)


Idents(Bos_ep_sub )<-Bos_ep_sub$Group
Bos_ep_sub.df <- FindMarkers(Bos_ep_sub, ident.1 = 'BOS', ident.2 = 'HC',recorrect_umi=FALSE)
Bos_ep_sub.df$gene <- rownames(Bos_ep_sub.df)

generate_Vol <- function(res.df){
  # New Volcano in ggplot
  up_genes <- res.df %>%
    dplyr::filter(gene %in% res.df$gene[res.df$avg_log2FC > 1.5 & res.df$p_val < 0.05])
  down_genes  <- res.df %>%
    dplyr::filter(gene %in% res.df$gene[res.df$avg_log2FC < -1.5 & res.df$p_val < 0.05])
  ns_genes_1 <- res.df %>%
    dplyr::filter(gene %in% res.df$gene[res.df$avg_log2FC < 1.5 &  res.df$avg_log2FC > -1.5])
  ns_genes_2 <- res.df %>%
    dplyr::filter(gene %in% res.df$gene[res.df$p_val > 0.05])
  # Legend
  res.df$Regulation <- NA
  res.df$Regulation[res.df$gene %in% up_genes$gene] <- 'UP'
  res.df$Regulation[res.df$gene %in% down_genes$gene] <- 'Down'
  res.df$Regulation[res.df$gene %in% ns_genes_1$gene] <- 'NS'
  res.df$Regulation[res.df$gene %in% ns_genes_2$gene] <- 'NS'
  plot<-ggplot(data=res.df, aes(x=avg_log2FC, y=-log10(p_val)))+
    geom_point(aes(colour = Regulation), 
               alpha = 0.2, 
               shape = 16,
               size = 2) + 
    geom_hline(yintercept=-log10(0.05), linetype='dashed')+
    geom_vline(xintercept = c(-1.5,1.5), linetype='dashed')+ 
    geom_point() +
    geom_point(data = up_genes,
               shape = 21,
               size = 4, 
               fill = "firebrick", 
               colour = "black") +
    geom_point(data = down_genes,
               shape = 21,
               size = 4, 
               fill = "blue", 
               colour = "black") +
    geom_point(data = ns_genes_1,
               shape = 21,
               size = 4, 
               fill = "green", 
               colour = "black") +
    geom_point(data = ns_genes_2,
               shape = 21,
               size = 4, 
               fill = "green", 
               colour = "black") +
    geom_label_repel(data = up_genes, # Add labels last to appear as the top layer
                     aes(label = gene), nudge_y = -0.4, nudge_x = 0.2, size=6, fontface='bold') +
    geom_label_repel(data = down_genes, # Add labels last to appear as the top layer
                     aes(label = gene), nudge_y = -0.4, nudge_x = 0.2, size=6, fontface='bold')+
    # geom_label_repel(data = ns_genes_1, # Add labels last to appear as the top layer  
    #                  aes(label = gene), nudge_y = -0.4, nudge_x = 0.2, size=6, fontface='bold')+
    scale_color_manual(values =  c("blue", "green",  "firebrick")) +
    theme_grey(base_size = 20)+
    theme(axis.text.x = element_text(face="bold"), text)
  
  print(plot)
}
setwd("~/Squidpy_Spatial_analysis/Airway/Airway_Mast_Basal/Volcano_plots_mast_basal_HC_vs_allconds")
pdf('basal_mast_COPDvHC_v2.pdf',15,12)
generate_Vol(CF_CXCL_ep_sub.df)
generate_Vol(Bos_ep_sub.df)
dev.off()

setwd("~/Squidpy_Spatial_analysis/acap_Mast_VolcanoPlots_all_layers_supl_figure")
pdf('acap_mast_COPDvHC.pdf',15,12)
generate_Vol(Mast_all_layers.df)
generate_Vol(acap_all_layers.df)
dev.off()

setwd("/media/haider/Files/Squidpy_Spatial_analysis/Mast_EP_VS_SUB")
pdf('mast_ep_vs_sub_ep.pdf', 10,10)
generate_Vol(mast_ep_sub.df)
dev.off()

setwd("/media/haider/Files/Squidpy_Spatial_analysis/Mast_EP_VS_SUB")
pdf('gata3_ep_vs_sub_ep.pdf', 10,10)
generate_Vol(gata3_sub_ep.df)
dev.off()


# CXCL Immune 
final_obj@active.ident <- factor(final_obj$Group)
# Dotplot to Determine that CXCL9 is most highly expressed in BOS

pdf('CXCL9_Dotplot_splitbyDiease.pdf', 5,5)
DotPlot(final_obj, features = 'CXCL9')
dev.off()

# Repeated the above procedure with a subseted group
immune_cxcl <- subset(final_obj, subset=subtype %in% c('cDC', 'Perivasc_macs', 'Mono/Macs', 'AM', 'CXCL+_Immune'))
immune_cxcl <- immune_cxcl  %>%
  RunPCA(npcs = 30, verbose = FALSE)
#ElbowPlot(immune_cxcl)
# one-liner to run Integration
Immune_cxcl_int = IntegrateLayers(object =immune_cxcl, method = "HarmonyIntegration", normalization.method = "SCT", verbose = T, new.reduction = "harmony", orig.reduction = "pca")
Immune_cxcl_int <- Immune_cxcl_int %>%  
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.2, verbose = FALSE)
saveRDS(Immune_cxcl_int,"~/Squidpy_Spatial_analysis/CXCL+_Immune/Immune_cxcl_int.rds")
FeaturePlot

# Annotation
new_IDs <- c('IM', 'AM', 'Mono','Perivasc_macs', 'CD68+', 'FCER1A+', 'CXCL_High_Lymph', 'cDC', 
             'pDC', 'CXCL_High', "IL1B+", 'CD163_High', 'CHIT1+',13,14,'AIF1+High', 'ALOX5AP_High', 'CD68_High',
             'Perivasc+FB',19, 'MS4A4A+', 21, 22, 23 )
# pdf('CXCL9_Dotplot_splitbyDiease_cdc_mono_am_perivasc_all_layers.pdf', 5,5)
# DotPlot(CDC_mono_perivasc_am, features = 'CXCL9')
# dev.off()
FeaturePlot(Immune_cxcl_int, cols=c('grey94', 'darkblue'), features='CXCL9', split.by = 'Group')
FeaturePlot(Immune_cxcl_int, cols=c('grey94', 'darkblue'), features='CXCL10', split.by = 'Group')
FeaturePlot(Immune_cxcl_int, cols=c('grey94', 'darkblue'), features='CXCL13', split.by = 'Group')

# Subclustering CXCL+_Immune cells
immune_cxcl_markers <- FindAllMarkers(Immune_cxcl_int, only.pos=TRUE, , recorrect_umi=FALSE)
Immune_cxcl_int@meta.data$subtype[Immune_cxcl_int@meta.data$seurat_clusters %in% c(9)] <- 'CXCL_high'
write.csv(immune_cxcl_markers, 'immune_cxcl_markers.csv')
# Repeat but for Airway Epitherlial and Sub-Epithelial
CDC_mono_perivasc_ep_am_sub <- subset(CDC_mono_perivasc_am, subset=Layers %in% c('Airway_Ep', 'Sub_Ep'))

pdf('CXCL9_Dotplot_splitbyDiease_cdc_mono_am_perivasc_ep_layers.pdf', 5,5)
DotPlot(CDC_mono_perivasc_ep_am_sub, features = 'CXCL9')
dev.off()

CDC_mono_perivasc_am@active.ident <- factor(CDC_mono_perivasc_am$ID)
CDC_mono_perivasc_am[["percent.cxcl9"]] <- PercentageFeatureSet(object = CDC_mono_perivasc_am, features = 'CXCL9')
data_frame_cxcl9  <- Percent_Expressing(CDC_mono_perivasc_am, features='CXCL9', group_by = 'ID')

Immune_cxcl_int@meta.data$subtype[Immune_cxcl_int@meta.data$seurat_clusters %in% c(9)]

CDC_mono_perivasc_am_avg <- AverageExpression(CDC_mono_perivasc_am, features = 'CXCL9')

# Boxplots
express_cxcl9  <- Percent_Expressing(CDC_mono_perivasc_am, features='CXCL9', group_by = 'ID')
data_frame_cxcl9 <- data.frame(ID=c('BOS11A','BOS12D','BOS13E','C48A','C50C','C51E','C52E','C56E','CF19','CF20','CF27', 'CF35','CF40','ID6','MTS12D','MTS21A','MTS27E','MTS29E','MTS48B','MTS7D'),
                                   Value=c(1.36594 ,1.332862, 5.736019, 1.337657, 2.978569, 2.169052, 3.354322, 0.8582463, 0.7769145,4.027043, 1.194997, 1.767141, 1.50021,
                                           1.416512, 2.871543, 0.3994469, 1.346336, 1.52263, 0.4859874, 1.913324))

data_frame_cxcl9[data_frame_cxcl9$ID %in% c('BOS11A', 'BOS12D', 'BOS13E'),]$ID <- 'BOS'
data_frame_cxcl9[data_frame_cxcl9$ID %in% c('MTS21A', 'MTS48B', 'MTS27E'),]$ID <- 'Asthma'
data_frame_cxcl9[data_frame_cxcl9$ID %in% c('ID6','MTS7D', 'MTS12D', 'MTS29E'),]$ID <- 'HC'
data_frame_cxcl9[data_frame_cxcl9$ID %in% c('C48A','C50C','C51E' ,'C52E', 'C56E'),]$ID <- 'COPD'
data_frame_cxcl9[data_frame_cxcl9$ID %in% c('CF19', 'CF20', 'CF27', 'CF35', 'CF40'),]$ID <- 'CF'

print(ggplot(data_frame_cxcl9, aes(x=ID, y=Value)) +
        geom_boxplot(outlier.shape = NA)+ geom_jitter(width = 0.2))


final_obj@active.ident <- final_obj$subtype
VlnPlot(final_obj, features='CXCL9', pt.size = 1)

# Annotation of CF34 and C53E
C48A <- readRDS("~/Downloads/COPD_gene_markers-selected/C48A.rds")
C51E <- readRDS("~/Downloads/COPD_gene_markers-selected/C51E.rds")
C53E <- readRDS("~/Downloads/COPD_gene_markers-selected/C53E.rds")

CF34 <- readRDS("~/Downloads/CF_gene_markers-selected/CF34.rds")
CF35 <- readRDS("~/Downloads/CF_gene_markers-selected/CF35.rds")
CF40 <- readRDS("~/Downloads/CF_gene_markers-selected/CF40.rds")

C53E@active.ident <- C53E$seurat_clusters
C53E <- RenameIdents(object = C53E,
                     '0'='FB',
                     '1'='gCap',
                     '2'='Gob',
                     '3'='Cil',
                     '4'='DES+_SMC',
                     '5'='CD4 Tcell',
                     '6'='Basal',
                     '7'='Sec',
                     '8'='LQ',
                     '9'= 'CD8 Tcell',
                     '10'='PLA2G2A+_FB',
                     '11'='APOD+_FB',
                     '12'='CDC',
                     '13'='GB',
                     '14' = 'Club',
                     '15'='LEC',
                     '16'='Mast',
                     '17'='SMC'
                     )

###########################Revised Object#######################################
# DotPlot
final_obj <- readRDS("/media/haider/Files/Squidpy_Spatial_analysis/final_annot_layers.rds")
# Remove low quality clusters 
final_obj <- subset(final_obj, subset=subtype %ni% c('LQ',  "AT2_EC_DB" ))
# Set markers

markers <- c('KRT7','SCEL','SFTPC','DMBT1','MUC5B','PIGR' ,'LTF','KRT5','KRT15', 'CD24', 'FOXJ1',
             'MMRN1', 'PROX1','KDR', 'CA4','ACKR1', 'PLVAP','RAMP2', 'FCN3','PDGFRB', 'KCNK3','ACTA2','MYH11', 'DES', 'LTBP2', 'FGFR4',
             'PLA2G2A','SFRP2', 'APOD','FBN1','CXCL13', 'CXCL10','CXCL9','KIT', 'MS4A2','S100A12','MMP9','MARCO','CD68',
             'CD14','MS4A4A', 'LYVE1','S100B','FGF2','CD19', 'MS4A1', 'FCMR','SEC11C',
             'MZB1', 'CD4','IL7R', 'GZMK', 'CD8A')
# Add Cxcl 9 high cell annotations
final_obj@active.ident <- factor(final_obj$subtype)
# cxcl_9_immune <- FeaturePlot(final_obj, cols=c('grey94', 'darkblue'), features = 'CXCL9')
cxcl9_selected_cells <- WhichCells(final_obj, expression = CXCL9 > 1.5)


final_obj$subtype <- as.character(final_obj$subtype)
final_obj$subtype[colnames(final_obj) %in% cxcl9_selected_cells] <- 'CXCL9_High'
final_obj$subtype <- final_obj@active.ident



# Fix T cell annotations
final_obj <- RenameIdents(final_obj, 'ThC'='CD4 Tcells')
final_obj <- RenameIdents(final_obj, 'TkC'='CD8 Tcells')

final_obj <- RenameIdents(final_obj, "CXCL+_Immune"='CXCL13_High')
final_obj$subtype <- final_obj@active.ident

levels(final_obj) <- c('CD8 Tcells', 'CD4 Tcells', 'GB', 'BC', 'cDC', 'Perivasc_macs', 'Mono/Macs','AM', 'PMN', 'Mast','CXCL9_High', "CXCL13_High",
                                                                  'APOD+_FB', 'PLA2G2A+_FB','FB','DES+_SMC','SMC', 'Peri',
                                                                  'FCN3+_EC', 'gCap','aCap','LEC', 
                                                                  'Cil', 'Basal', 'Sec', 'Gob','DMBT1+_AT2', 'AT2', 'AT1')

# Create H5ad object
DefaultAssay(final_obj ) <- "Xenium"

final_obj  <- JoinLayers(
  final_obj ,
  assay = "Xenium"
)

fovs <- Images(final_obj)

coords_list <- lapply(fovs, function(fov) {
  
  tmp <- GetTissueCoordinates(
    final_obj,
    image = fov
  )
  
  tmp$fov <- fov
  
  tmp
})

coords <- do.call(rbind, coords_list)
# -----------------------------
# 1. Convert Seurat -> SCE
# -----------------------------
sce <- as.SingleCellExperiment(
  final_obj,
  assay = "Xenium"
)
# Match coordinates to SCE cell order
idx <- match(colnames(sce), coords$cell)

# Create spatial coordinate matrix
spatial <- as.matrix(coords[idx, c("x", "y")])

# Make sure rownames match the cells
rownames(spatial) <- colnames(sce)

# Add to reducedDims
reducedDim(sce, "spatial") <- spatial

writeH5AD(
  sce,
  file = "~/Whole_Tissue_All_Layers_Lung.h5ad",
  X_name = 'counts'
)
# -----------------------------
# 5. Add existing reductions
# -----------------------------
reducedDim(sce, "X_pca") <- Embeddings(final_obj, "pca")
reducedDim(sce, "X_harmony") <- Embeddings(final_obj, "harmony")
reducedDim(sce, "X_umap") <- Embeddings(final_obj, "umap")

# -----------------------------
# 6. Write AnnData
# -----------------------------
writeH5AD(
  sce,
  file = "xenium_nico.h5ad",
  X_name = "counts"
)

################################################################################
setwd("/media/haider/Files/Squidpy_Spatial_analysis/ALL_Sample_DOTPLOT")
library(circlize)
cluster_rand <- rand_color(length(levels(factor(final_obj$subtype))))
set.seed(1)

final_obj$subtype <- final_obj@active.ident

pdf('all_sample_dotplot_and_dimplot_new.pdf',20, 10)
DotPlot(final_obj, features = markers,  cols = c("blue", "red")) + RotatedAxis()
DimPlot(final_obj, group.by = 'subtype', label=T)
dev.off()


pdf('all_sample_dotplot_new_v2.pdf',20, 10)
DotPlot(final_obj, features = markers,  cols = c("blue", "red")) + RotatedAxis()
dev.off()

pdf('all_sample_dimplot_new.pdf',20, 10)
DimPlot(final_obj, group.by = 'subtype', label=T, cols = cluster_rand, raster = FALSE)
dev.off()

### Whole Tissue Plot1cell Dotplots for cell


plotcell_gene_markers <- function(obj, cell_type) {
  obj <- subset(obj, subset = subtype %in% c(cell_type))
  markers <- FindAllMarkers(obj, only.pos = TRUE, recorrect_umi=FALSE)
  write.csv(markers, paste0('/media/haider/Files/Squidpy_Spatial_analysis/Plot1cell_Dotplot_whole_tissue/','Mono', '_markers.csv'))
}

final_obj@active.ident <- factor(final_obj$Group)
final_obj <- RenameIdents(final_obj, 'BOS'='CLAD')
levels(final_obj)<-c('HC', 'COPD', 'Asthma', 'CF', 'CLAD')
final_obj$Group <- final_obj@active.ident
final_obj@active.ident <- factor(final_obj$subtype)

plotcell_gene_markers(final_obj, 'Basal')
plotcell_gene_markers(final_obj, 'Sec')
plotcell_gene_markers(final_obj, 'GB')
plotcell_gene_markers(final_obj, 'Mast')
plotcell_gene_markers(final_obj, 'Mono/Macs')
plotcell_gene_markers(final_obj, 'PLA2G2A+_FB')
plotcell_gene_markers(final_obj, 'FB')
plotcell_gene_markers(final_obj, 'gCap')
plotcell_gene_markers(final_obj, 'aCap')
plotcell_gene_markers(final_obj, 'BC')
plotcell_gene_markers(final_obj, 'LEC')

library(ggplot2)
library(Seurat)
library(reshape2)
library(scales)

# Helper function to extract counts safely in Seurat v5
# Replaces old internal extraction logic with standard FetchData
.extract_v5_counts <- function(seu_obj, features, cell.types, meta.groups, splitby = NULL) {
  # Standardize metadata pull
  vars_to_fetch <- c(features, meta.groups)
  if (!is.null(splitby)) {
    vars_to_fetch <- unique(c(vars_to_fetch, splitby))
  }
  
  # Fetch data uses the active assay's 'data' layer automatically in v5
  data_mat <- FetchData(object = seu_obj, vars = vars_to_fetch, clean = FALSE)
  data_mat$celltype <- Idents(seu_obj)
  
  # Filter for requested cell types
  data_mat <- data_mat[data_mat$celltype %in% cell.types, , drop = FALSE]
  return(data_mat)
}

complex_dotplot_single <- function(
    seu_obj, 
    feature, 
    celltypes=NULL,
    groups,
    splitby=NULL,
    color.palette = NULL,
    font.size = 12,
    strip.color=NULL,
    do.scale=TRUE,
    scale.by='radius',
    vmin = NULL, 
    vmax = NULL
){
  .clip_vals <- function(x, vmin = NULL, vmax = NULL){
    if(!is.null(vmin)) x[x < vmin] <- vmin
    if(!is.null(vmax)) x[x > vmax] <- vmax
    x
  }
  
  if(is.null(color.palette)){
    color.palette <- colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255)
  }
  
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  
  # Seurat v5 compatibility update for fetch identity levels
  if(is.null(celltypes)){
    celltypes <- levels(Idents(seu_obj))
  }
  
  if(length(groups) == 1){
    groups_level <- levels(factor(seu_obj@meta.data[, groups]))
    if (is.null(groups_level)){
      seu_obj@meta.data[, groups] <- factor(seu_obj@meta.data[, groups], levels = names(table(seu_obj@meta.data[, groups])))
      groups_level <- levels(seu_obj@meta.data[, groups])
    } 
    
    if(!is.null(splitby)){
      if (is.null(levels(factor(seu_obj@meta.data[, splitby])))){
        seu_obj@meta.data[, splitby] <- factor(seu_obj@meta.data[, splitby], levels = names(table(seu_obj@meta.data[, splitby])))
      }
      splitby_level <- levels(factor(seu_obj@meta.data[, splitby]))
      
      # V5 Adapted extraction
      count_df <- .extract_v5_counts(seu_obj, features = feature, cell.types = celltypes, meta.groups = groups, splitby = splitby)
      
      count_df$new_group <- paste(count_df[, groups], count_df[, "celltype"], count_df[, splitby], sep = "___")
      exp_df <- aggregate(as.formula(paste(feature, "~ new_group")), data=count_df, FUN=function(x){mean(expm1(x))}) 
      pct_df <- aggregate(as.formula(paste(feature, "~ new_group")), data=count_df, FUN=function(x){length(x[x > 0]) / length(x)})
      
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      data_plot <- merge(exp_df, pct_df, by='new_group')
      
      data_plot$groups <- sapply(strsplit(data_plot$new_group, split = "___"), function(x) x[1])
      data_plot$celltype <- sapply(strsplit(data_plot$new_group, split = "___"), function(x) x[2])
      data_plot$splitby <- sapply(strsplit(data_plot$new_group, split = "___"), function(x) x[3])
      
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$splitby <- factor(data_plot$splitby, levels = splitby_level)
      data_plot$celltype <- factor(data_plot$celltype, levels = rev(celltypes))
    } else {
      # V5 Adapted extraction
      count_df <- .extract_v5_counts(seu_obj, features = feature, cell.types = celltypes, meta.groups = groups)
      
      count_df$new_group <- paste(count_df[, groups], count_df[, "celltype"], sep = "___")
      exp_df <- aggregate(as.formula(paste(feature, "~ new_group")), data=count_df, FUN=function(x){mean(expm1(x))})
      pct_df <- aggregate(as.formula(paste(feature, "~ new_group")), data=count_df, FUN=function(x){length(x[x > 0]) / length(x)})
      
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      data_plot <- merge(exp_df, pct_df, by='new_group')
      
      data_plot$groups <- sapply(strsplit(data_plot$new_group, split = "___"), function(x) x[1])
      data_plot$celltype <- sapply(strsplit(data_plot$new_group, split = "___"), function(x) x[2])
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$celltype <- factor(data_plot$celltype, levels = rev(celltypes))
    }
    
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 2)
    data_plot$avg.exp <- as.vector(scale(data_plot$avg.exp))
    data_plot$avg.exp <- .clip_vals(data_plot$avg.exp, vmin, vmax)
    
    p <- ggplot(data_plot, aes(y = celltype, x = groups)) +  
      geom_tile(fill="white", color="white") +
      geom_point(aes(colour=avg.exp, size =pct.exp))  +  
      scale_color_gradientn(colours  =  color.palette )+ 
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = font.size),
            plot.title = element_text(size = (font.size +2), hjust = 0.5, face = 'bold'),
            axis.text = element_text(size = font.size),
            legend.text=element_text(size=(font.size-2)),
            legend.key = element_rect(fill = NA, colour = NA) ,
            legend.title = element_text(size = (font.size)),
            strip.text = element_text( size = font.size),
            legend.position="right")+
      ylab("")+xlab("")+ggtitle(feature)+ guides(
        size = guide_legend(override.aes = list(shape = 21, fill = NA, colour = "black"))
      )
    
    if(do.scale){
      p = p + scale_size(range = c(0, 10))
    } else {
      if(max(data_plot$pct.exp) >= 20){
        p = p + scale_size(range = c(0, 10))
      } else {
        p = p + scale.func(range = c(0, 10), limits = c(0, 20))
      }
    }
    
    if(!is.null(splitby)){
      p <- p + facet_wrap(~splitby, scales = 'free_x')
      g <- change_strip_background(p, type = 'top', strip.color = strip.color)
      grid::grid.draw(g)
    } else {
      return(p)
    }
    
  } else { ### Multiple groups check
    gene_count <- .extract_v5_counts(seu_obj=seu_obj, features = feature, cell.types = celltypes, meta.groups = groups, splitby = splitby)
    allgroups <- c(groups, splitby)
    
    for(i in 1:length(allgroups)){
      if (is.null(levels(factor(seu_obj@meta.data[, allgroups[i]])))){
        seu_obj@meta.data[, allgroups[i]] <- factor(seu_obj@meta.data[, allgroups[i]], levels = names(table(seu_obj@meta.data[, allgroups[i]])))
      }
      group_level <- levels(factor(seu_obj@meta.data[, allgroups[i]]))
      gene_count[, allgroups[i]] <- factor(gene_count[, allgroups[i]], levels = group_level)
    }
    
    gene_count$celltype <- factor(gene_count$celltype, levels = celltypes)
    all_levels <- list()
    for(i in 1:length(groups)){
      if (is.null(levels(factor(seu_obj@meta.data[, groups[i]])))){
        seu_obj@meta.data[, groups[i]] <- factor(seu_obj@meta.data[, groups[i]], levels = names(table(seu_obj@meta.data[, groups[i]])))
      }
      group_level <- levels(factor(seu_obj@meta.data[, groups[i]]))
      all_levels[[i]] <- group_level
    }
    all_levels <- as.character(unlist(all_levels))
    
    data_plot <- list()
    for(i in 1:length(groups)){
      count_df <- gene_count
      count_df$new_group <- paste(gene_count[, groups[i]], gene_count[, "celltype"], sep = "___")
      
      exp_df <- aggregate(as.formula(paste(feature, "~ new_group")), data=count_df, FUN=function(x){mean(expm1(x))})
      pct_df <- aggregate(as.formula(paste(feature, "~ new_group")), data=count_df, FUN=function(x){length(x[x > 0]) / length(x)})
      
      colnames(exp_df)[2] <- "avg.exp"
      colnames(pct_df)[2] <- "pct.exp"
      df1 <- merge(exp_df, pct_df, by='new_group')
      df1$groupID <- groups[i]
      data_plot[[i]] <- df1
    }
    data_plot <- do.call("rbind", data_plot)
    
    data_plot$groups <- sapply(strsplit(data_plot$new_group, split = "___"), function(x) x[1])
    data_plot$celltype <- sapply(strsplit(data_plot$new_group, split = "___"), function(x) x[2])
    data_plot$groups <- factor(data_plot$groups, levels = all_levels)
    data_plot$celltype <- factor(data_plot$celltype, levels = rev(celltypes))
    data_plot$groupID <- factor(data_plot$groupID, levels = groups)
    
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 2)
    data_plot$avg.exp <- as.vector(scale(data_plot$avg.exp))
    data_plot$avg.exp <- .clip_vals(data_plot$avg.exp, vmin, vmax)
    
    if(is.null(splitby)){
      p <- ggplot(data_plot, aes(y = celltype, x = groups)) +  
        geom_tile(fill="white", color="white") +
        geom_point(aes(colour=avg.exp, size =pct.exp))  +  
        scale_color_gradientn(colours  =  color.palette )+ 
        theme(panel.background = element_rect(fill = "white", colour = "black"),
              axis.text.x = element_text(angle = 45, hjust = 1, size = font.size),
              plot.title = element_text(size = (font.size +2), hjust = 0.5, face = 'bold'),
              axis.text = element_text(size = font.size),
              legend.text=element_text(size=(font.size-2)),
              legend.title = element_text(size = (font.size)),
              legend.key = element_rect(fill = NA, colour = NA) ,
              strip.text = element_text( size = font.size),
              legend.position="right")+
        ylab("")+xlab("")+ggtitle(feature)+facet_wrap(~groupID, scales = 'free_x')+ guides(
          size = guide_legend(override.aes = list(shape = 21, fill = NA, colour = "black"))
        )
      if(do.scale){
        p = p + scale_size(range = c(0, 10))
      } else {
        if(max(data_plot$pct.exp) >= 20){
          p = p + scale_size(range = c(0, 10))
        } else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 20))
        }
      }
      g <- change_strip_background(p, type = 'top',  strip.color = strip.color)
      grid::grid.draw(g)
    } else {
      df2 <- reshape2::melt(gene_count[, c(groups, splitby)], id.vars = splitby, measure.vars = groups)
      df2 <- df2[!duplicated(df2$value),]
      colnames(df2)[colnames(df2) == "value"] <- "groups"
      
      data_plot2 <- list()
      for(i in 1:length(groups)){
        df3 <- data_plot[data_plot$groupID == groups[i],]
        df4 <- df2[df2$variable == groups[i], c('groups', splitby[i])]
        colnames(df4)[2] <- "split"
        df5 <- merge(df3, df4, by='groups')
        data_plot2[[i]] <- df5
      }
      data_plot2 <- do.call("rbind", data_plot2)
      
      fill_x1 <- grDevices::rainbow(length(groups), alpha = 0.5)
      fill_x2 <- list()
      for(i in 1:length(splitby)){
        n_col <- unique(gene_count[, splitby[i]])
        fill_x2[[i]] <- scales::hue_pal(l=90)(length(n_col))
      }
      fill_x2 <- as.character(unlist(fill_x2))
      fill_x <- c(fill_x1, fill_x2)
      
      p <- ggplot(data_plot2, aes(y = celltype, x = groups)) +  
        geom_tile(fill="white", color="white") +
        geom_point(aes(colour=avg.exp, size =pct.exp))  +  
        scale_color_gradientn(colours  =  color.palette )+ 
        theme(panel.background = element_rect(fill = "white", colour = "black"),
              axis.text.x = element_text(angle = 45, hjust = 1, size = font.size),
              plot.title = element_text(size = (font.size +2), hjust = 0.5, face = 'bold'),
              axis.text = element_text(size = font.size),
              legend.text=element_text(size=(font.size-2)),
              legend.key = element_rect(fill = NA, colour = NA) ,
              legend.title = element_text(size = (font.size)),
              strip.text = element_text( size = font.size),
              legend.position="right")+
        ylab("")+xlab("")+ggtitle(feature)+
        ggh4x::facet_nested(~ groupID + split, scales = "free_x", 
                            strip = ggh4x::strip_nested(background_x = ggh4x::elem_list_rect(fill = fill_x)))+ 
        guides(size = guide_legend(override.aes = list(shape = 21, fill = NA, colour = "black")))
      
      if(do.scale){
        p = p + scale_size(range = c(0, 10))
      } else {
        if(max(data_plot$pct.exp) >= 20){
          p = p + scale_size(range = c(0, 10))
        } else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 20))
        }
      }
      return(p)
    }
  }
}

complex_dotplot_multiple <- function(
    seu_obj, 
    features, 
    celltypes=NULL,
    groups, 
    color.palette = NULL,
    strip.color = NULL,
    vmin = NULL, 
    vmax = NULL,
    strip.alpha = 1
){
  # Added a simple fallback loop handler if progress_bar library isn't active
  if(requireNamespace("progress", quietly = TRUE)) {
    pb <- progress::progress_bar$new(
      format = "  Plotting [:bar] :percent eta: :eta",
      clear = FALSE, total = length(features), width = 100)
  } else {
    pb <- NULL
  }
  
  plot_list <- list()
  for(i in 1:length(features)){
    # Prevent multi-threading conflicts during individual V5 matrix calls
    pp <- tryCatch({
      complex_dotplot_single(seu_obj = seu_obj, feature = features[i], groups = groups, celltypes = celltypes, vmin = vmin, vmax = vmax)
    }, error = function(e) { return(NULL) })
    
    if(!is.null(pp)) {
      # Extract underlying structure
      plot_data <- pp
      if(inherits(pp, "ggplot")) { plot_data <- pp$data }
      plot_data$gene <- features[i]
      plot_list[[i]] <- plot_data
    }
    
    if(!is.null(pb)) pb$tick()
    Sys.sleep(0.01)
  }
  
  all_data <- do.call('rbind', plot_list)
  all_data$gene <- factor(all_data$gene, levels = rev(features)) 
  
  if(is.null(celltypes)) {
    celltypes <- levels(Idents(seu_obj))
  }
  all_data$celltype <- factor(all_data$celltype, levels = celltypes)
  
  if(is.null(color.palette)){
    color.palette <- colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255)
  }
  
  p <- ggplot(all_data, aes(x = groups, y = gene)) +   
    geom_tile(fill="white", color="white") +
    geom_point(aes(colour=avg.exp, size =pct.exp), alpha=0.9)  +   
    scale_color_gradientn(colours  =  color.palette)+  
    scale_size(range = c(0, 10))+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size=16),
          plot.title = element_text(size = 16,hjust = 0.5, face = 'bold'),
          axis.text = element_text(size = 12),
          axis.title=element_text(size=8),
          legend.text=element_text(size=8),
          legend.key = element_rect(fill = NA, colour = NA) ,
          legend.title = element_text(size = 12),
          legend.position="right",
          strip.text = element_text(size = 14,colour = 'black',face = 'bold'))+
    ylab("")+xlab("")+ggtitle('')+
    facet_wrap(~celltype, ncol = length(celltypes)) + 
    guides(size = guide_legend(override.aes = list(shape = 21, fill = NA, colour = "black")))
  
  if (!is.null(strip.color)) strip.color <- scales::alpha(strip.color, strip.alpha)
  g <- change_strip_background(p, type = 'top', strip.color = strip.color)
  grid::grid.draw(g)
}

change_strip_background <- function(
    ggplt_obj, 
    type = "top",
    strip.color=NULL
){
  g <- ggplot_gtable(ggplot_build(ggplt_obj))
  if(type == "top"){
    strip_both <- which(grepl('strip-t', g$layout$name))
    fills<-strip.color
    if(is.null(fills)){
      fills<- scales::hue_pal(l=90)(length(strip_both))
    }
  } else if(type=="right"){
    strip_both <- which(grepl('strip-r', g$layout$name))
    fills<-strip.color
    if(is.null(fills)){
      fills<- scales::hue_pal(l=90)(length(strip_both))
    }
  } else {
    strip_t <- which(grepl('strip-t', g$layout$name))
    strip_r <- which(grepl('strip-r', g$layout$name))
    strip_both<-c(strip_t, strip_r)
    fills<-strip.color
    if(is.null(fills)){
      fills<- c(scales::hue_pal(l=90)(length(strip_t)),scales::hue_pal(l=90)(length(strip_r)))
    }
  }
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  g
}

plot1cell_dotplot1 <- function(obj, cell_type, markers, color.palette = NULL) {
  
  meta <- obj@meta.data
  
  stopifnot("subtype" %in% colnames(meta))
  
  # SAFE SUBSET (NO WhichCells, NO Idents)
  cells.use <- rownames(meta)[meta$subtype %in% cell_type]
  
  if (length(cells.use) == 0) {
    stop("No cells found for subtype")
  }
  
  object <- obj[, cells.use]
  
  if (is.null(color.palette)) {
    color.palette <- colorRampPalette(
      c('#0200ad', '#fbfcbd', '#ff0000')
    )(255)
  }
  
  if (cell_type == 'Mono/Macs'){
    cell_type <- 'Mono_Macs'
  }
  pdf( paste0(cell_type, "_Dot_Disease.pdf"), 5,8)
  
  plots <- complex_dotplot_multiple(
    object,
    features = markers,
    groups = "Group",
    color.palette = color.palette
  )
  
  for (p in plots) print(p)
  
  dev.off()
}

final_obj@active.ident <- final_obj$subtype
acap_markers <-  c("S100A4",'TMEM100', 'CA4', 'PECAM1','SHANK3',
                   'GNG11', 'MYO6', 'ACE', 'CD34', 'KDR', 
                    'APOLD1', 'MTUS1', 'BAIAP2L1', 'PDCD1LG2',
                   'F3', 'SMPD3', 'IGFBP3')
plot1cell_dotplot1(final_obj, 'aCap', acap_markers)

LEC_markers <-  c("MMRN1", "TM4SF18", "GNG11", "KDR", "S100A4","PDPN",
                   "CCN2", "MTUS1", "STC1", 'FAM184A', 
                  'BMX', 'NCEH1', 'PI3',
                   'THY1', 'SMIM24')
plot1cell_dotplot1(final_obj, 'LEC', LEC_markers)

gCap_markers <-  c("ECSCR", "SNX3","S100A4", "ITGB4", "CD34", "MALL", 
                   "CCN2", "NTRK21", "PCOLCE2", "IL1A",
                    "PROX1","APOLD1", "MTUS1", "MYC", "PLVAP1", "STC1",
                    'SLC7A11')
plot1cell_dotplot1(final_obj, 'gCap', gCap_markers)

BC_markers <-  c('BANK1', "MS4A1", 'FCMR', "IRF8", "CD79A", "CD1C",
                 "GPR183", "GNG11", "MYC", "TRPC6", "MPEG1", 
                 "IL33A1", "CXCL13",
                 "AQP9", "GZMB", "LAG3", "GPR171", "IL4")
plot1cell_dotplot1(final_obj, 'BC', BC_markers)

GB_markers <-  c("CD63", "S100A4", "SEC11C", "MZB1", "PIM2",
                 "ADAM28", "P2RX1","MYC", "CYSLTR1",
                 "FCMR", "CD79A", "TNFRSF13C")
plot1cell_dotplot1(final_obj, 'GB', GB_markers)

mast_markers <- c("ALOX5", "PTGS1", "HPGDS", "IL1RL1",
                  "FCER1A", "GPR183", "IL5","CD38", "IL4",
                  "MKI67", "IL13", "AREG")
plot1cell_dotplot1(final_obj, 'Mast', mast_markers)


sec_markers <- c("KRT7", "PI3", "SERPINA3", "LTF",
                 'EHF', "KLF5", "KDR", "CSTA", "RARRES1","CFB",
                 "CXCL6", "ELF3", 'MUC5B', 'BPIFA1', 'MUC1', "MALL", "PEBP4",
                 "DUOX1", "LTA4H")
plot1cell_dotplot1(final_obj, 'Sec', sec_markers)

fb_markers <-  c("SERPINA3", "VIM", "S100A4", "FGF7", "HIF1A", "CD63",
                 "CCN2", "SLIT3", "FOS", "SNAI2", "ANO1", "PCOLCE2", 'APOD',
                 "IGFBP3", "SFRP2", "IGF1", "STEAP4", "TGFB3",
                 "THY1", "RND1", "PIM1")
plot1cell_dotplot1(final_obj, 'FB', fb_markers)

basal_markers <- c("KLF5", "EGFR", "IL33", "FOS", "TP63", "GSDMC", "KRT5",
                   "ANO1", "SCEL", "MKI67",'KRT14', "CXCL6",
                   'PIGR', 'SOX9', "THBS2", "ADAMTS1"
)
plot1cell_dotplot1(final_obj, 'Basal', basal_markers)

PLA2G2A_FB_markers <- c("SERPINA3", "FGF7", "S100A4", "PLA2G2A", "VIM",
                        "PCOLCE2", "CCN2", "MFAP5", 'CD34', "SNAI2", "GSDME",
                        "APOD","SFRP2", "THY1", "THBS2",
                        "LTBP2", 'COL5A2', 'ACTA2')
plot1cell_dotplot1(final_obj, 'PLA2G2A+_FB', PLA2G2A_FB_markers )

mono_markers <- c("S100A4", "CD14", "CD68", "CD300E", "MPEG1",
                  "CD1C", "CLEC10A","CD1A", "CHIT1",
                  "CLEC4E", "ADAM28", "LYVE1",
                   "CXCL9"
                  )
plot1cell_dotplot1(final_obj, 'Mono/Macs', mono_markers)

## Dotplot for genes per condition
final_obj@active.ident <- final_obj$Group
pdf('IL33_COPD_Dotplot', 5,5)
DotPlot(final_obj, features = 'IL33')
dev.off()


pdf('CXCL13_Dotplot', 5,5)
DotPlot(final_obj, features = 'CXCL13')
dev.off()

pdf('IL5_Dotplot', 5,5)
DotPlot(final_obj, features = 'IL5')
dev.off()

##convert seurat to h5ad object for R-L analysis

library(SingleCellExperiment)
library(zellkonverter)

obj <- final_obj
DefaultAssay(obj) <- "Xenium"

obj[["SCT"]] <- NULL
obj[["BlankCodeword"]] <- NULL
obj[["ControlCodeword"]] <- NULL
obj[["ControlProbe"]] <- NULL

obj <- JoinLayers(obj)
sce <- as.SingleCellExperiment(obj)

writeH5AD(
  sce,
  "/media/haider/Files/Squidpy_Spatial_analysis/annot_v3_xenium.h5ad"
)

# Featureplots for Plasma cells
pdf("~/Squidpy_Spatial_analysis/Airway/Featureplots_Sept_2026/whole_tissue_split_disease_CD27.pdf",10, 4)
FeaturePlot(final_obj, raster = T, features = 'CD27', split.by = 'Group', assay = 'SCT')
dev.off()
pdf("~/Squidpy_Spatial_analysis/Airway/Featureplots_Sept_2026/whole_tissue_split_disease_CD38.pdf",10, 4)
FeaturePlot(final_obj, raster = T, features = 'CD38', split.by = 'Group', assay = 'SCT')
dev.off()

pdf("~/Squidpy_Spatial_analysis/Airway/Featureplots_Sept_2026/whole_tissue_split_disease_MZB1.pdf",10, 4)
FeaturePlot(final_obj, raster = T, features = 'MZB1', split.by = 'Group', assay = 'SCT')
dev.off()

pdf("~/Squidpy_Spatial_analysis/Airway/Featureplots_Sept_2026/whole_tissue_split_disease_TNFRSF17.pdf",10, 4)
FeaturePlot(final_obj, raster = T, features = 'TNFRSF17', split.by = 'Group', assay = 'SCT')
dev.off()

pdf("~/Squidpy_Spatial_analysis/Airway/Featureplots_Sept_2026/whole_tissue_split_disease_TNFRSF13B.pdf",10, 4)
FeaturePlot(final_obj, raster = T, features = 'TNFRSF13B', split.by = 'Group', assay = 'SCT')
dev.off()

# Save complete object
saveRDS(final_obj, "~/Squidpy_Spatial_analysis/spatial_obj_whole_lung.rds")



