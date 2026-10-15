library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
library(anndata)
library(plot1cell)
library(SeuratExtend)
library(ggrepel)
library(anndata)
library(scCustomize)
plan(multisession(workers = 20))
options(future.globals.maxSize = 8000 * 1024^3)
color.palette = colorRampPalette(c('#0200ad','#fbfcbd','#ff0000'))(255)
setwd("~/Spatial_Data")
set.seed(12)
################################################################################
data <- read_h5ad("~/Squidpy_Spatial_analysis/Airway_subclass_obj")
`%ni%` <- Negate(`%in%`)

final_obj <- readRDS("/media/haider/Files/Spatial_Data/annotated_comb_final.rds")
# Remove low quality clusters 
final_obj <- subset(final_obj, subset=subtype %ni% c('LQ',  "AT2_EC_DB" ))
# Get genelists for all genes
all_genes <- final_obj %>% PrepSCTFindMarkers() %>% FindAllMarkers(only.pos = TRUE)

#FeaturePlots
copd_hc_subset <- subset(final_obj, Group %in% c('HC', 'COPD'))

pdf('COPD_HC_DimPlot.pdf', 10,10)
DimPlot(copd_hc_subset, label=T, order=T, raster = F, cols = cluster_rand)
dev.off()
pdf('COPD_HC_IL1RL1_IL33_Featureplot.pdf',8,8)
FeaturePlot(copd_hc_subset,features = c('IL1RL1', 'IL33'), split.by = 'Group', raster = F)
dev.off()

# Dimplot for all Disease Groups
cluster_rand <- rand_color(length(levels(factor(final_obj$Group))))
setwd("~/Squidpy_Spatial_analysis/DimPlot_ALL_DIseaseGroups_Combined_Data_copd_hc_bos_asthma_cf")
pdf('All_Disease_Groups_DimPlot.pdf', 8,8)

DimPlot(final_obj_down, raster = F, group.by = 'Group', shuffle = TRUE,  cols = c('orange', 'darkblue', 'pink', 'red', 'cyan'), sizes.highlight = 4)
dev.off()


setwd("~/Squidpy_Spatial_analysis/DimPlot_ALL_DIseaseGroups_Combined_Data_copd_hc_bos_asthma_cf")

final_obj_down <- subset(final_obj, downsample=10000)
cluster_rand <- rand_color(length(levels(final_obj_down)))
pdf('Allgroups_DimPlot.pdf', 10,10)
DimPlot(final_obj_down, raster = F, order=T, cols = cluster_rand, label=T, repel = T, label.size = 3, label.color = 'black')
dev.off()

#Layers
ep <- data$obs[data$obs$Layers %in% 'Airway_Ep',]
sub <- data$obs[data$obs$Layers %in% 'Sub_Ep',]
smg <- data$obs[data$obs$Layers %in% 'SMG',]

final_obj$X <- sapply(strsplit(as.character(colnames(final_obj)), "_"), `[`, 2)
ep_obj <- subset(final_obj, subset = X %in% ep$cell_ids)
ep_obj@active.assay <- 'SCT'
ep_obj <- subset(ep_obj, subset= ID %ni% c('MTS43D', 'MTS38D'))


sub_obj <- subset(final_obj, subset = X %in% sub$cell_ids)
sub_obj@active.assay <- 'SCT'
sub_obj <- subset(sub_obj, subset= ID %ni% c('MTS43D', 'MTS38D'))


smg_obj <- subset(final_obj, subset = X %in% smg$cell_ids)
smg_obj@active.assay <- 'SCT'
smg_obj <- subset(smg_obj, subset= ID %ni% c('MTS43D', 'MTS38D'))


basal <- subset(airway_obj, subset = subtype %in% c("Basal"))

# Generating Dotplots for Genes in the Basal subset across Disease Conditions
DGE_list <- c('FSCN1','SOX9','FOS','THBS2','ADAMTS1','S100B','CCR7','CXCL9','CXCL6', 'CFB', 'CXCL5', 'S100A7', 'TP73', 'ANO1', 'MKI67', 'CDK1', 'TOP2A', 'FOXJ1')

DGE_list2 <- c('IL33FULL-exon7-exon8', 'LTBP2','IL33A1-exon1A1-exon2', 'IL33')

pdf('Allgenes_basal.pdf', 20, 20)
print(DotPlot(basal, features = DGE_list2, split.by = 'Group', cols = c('red', 'black', 'darkgreen', 'blue', 'brown')))

print(DotPlot(basal, features = DGE_list, split.by = 'Group', cols = c('red', 'black', 'darkgreen', 'blue', 'brown')))
dev.off()

# BoxPlots

library(ggplot2)
library(ggrepel)
ep_obj <- readRDS("~/Squidpy_Spatial_analysis/ep_obj.rds")
create_table <- function(obj){
  final_table2 <- table(obj@active.ident, obj$ID)
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
  #final_table2[,20] <- final_table2[,20]/sum(final_table2[,20])*100
  final_table2 <- as.data.frame(final_table2)
  final_table2$Var2 <- as.character(final_table2$Var2)
  final_table2[final_table2$Var2 %in% c('BOS11A', 'BOS12D', 'BOS13E'),]$Var2 <- 'BOS'
  final_table2[final_table2$Var2 %in% c('MTS21A', 'MTS48B', 'MTS27E'),]$Var2 <- 'Asthma'
  final_table2[final_table2$Var2 %in% c('ID6','MTS7D', 'MTS12D', 'MTS29E', 'MTS38D', 'MTS43D'),]$Var2 <- 'HC'
  final_table2[final_table2$Var2 %in% c('C48A','C50C','C51E' ,'C52E', 'C56E'),]$Var2 <- 'COPD'
  final_table2[final_table2$Var2 %in% c('CF19', 'CF20', 'CF27', 'CF35', 'CF40'),]$Var2 <- 'CF'
  final_table2$Var2 <- as.factor(final_table2$Var2)
  final_table2$Var2 <- factor(final_table2$Var2, levels = c('HC', 'Asthma', 'BOS', 'CF', 'COPD'))
  return(final_table2)
  }
ep_table <- create_table(ep_obj)
sub_table <- create_table(sub_obj)
smg_table <- create_table(smg_obj)

write.csv(smg_table, 'smg_freq_table.csv')
write.csv(sub_table, 'sub_freq_table.csv')

immune_list <- c('TkC', 'ThC', 'cDC', 'AM', 'Mast', 'GB', 'BC', 'PMN',"Perivasc_macs", "CXCL+_Immune", 'Mono/Macs')
alveolar_list <- c('AT1', 'AT2', 'AT2_EC_DB', 'DMBT1+_AT2', 'aCap', 'gCap','FCN3+_EC')
vascular_list <- c('Basal', 'Cil', 'Sec', 'Gob', 'SMC', 'DES+_SMC', 'PLA2G2A+_FB')
mes_list <- c('APOD+_FB','FB','Peri', 'LEC')
all_cells <- append(immune_list, c(alveolar_list, vascular_list, mes_list))

immune_freq_table <- ep_table[ep_table$Var1 %in% immune_list,]
vascular_freq_table <- ep_table[ep_table$Var1 %in% vascular_list,]



write.csv(immune_freq_table, "~/Squidpy_Spatial_analysis/Airway/Airway_EP_freq/immune.csv")
write.csv(vascular_freq_table, "~/Squidpy_Spatial_analysis/Airway/Airway_EP_freq/vascular.csv")


whisker <- function(sample_list, obj){
  for (i in 1:length(sample_list)){
    print(ggplot(obj[obj$Var1 %in% c(sample_list[i]),], aes(x=`Var2`, y=Freq)) +
            geom_boxplot(outlier.shape = NA)+ geom_jitter(width = 0.2)+ggtitle(sample_list[i]))
  }}
setwd("~/Squidpy_Spatial_analysis/Airway/BoxPlots")
# Epithelial
pdf('immune_ep_reorder.pdf',8,8)
whisker(immune_list, ep_table)
dev.off()
pdf('vascular_ep.pdf',8,8)
whisker(vascular_list, ep_table)
dev.off()
pdf('mes_ep.pdf',8,8)
whisker(mes_list, ep_table)
dev.off()

# Sub-Epithelial
pdf('immune_sub.pdf',8,8)
whisker(immune_list, sub_table)
dev.off()
pdf('vascular_sub.pdf',8,8)
whisker(vascular_list, sub_table)
dev.off()
pdf('mes_sub.pdf',8,8)
whisker(mes_list, sub_table)
dev.off()

# SMG
pdf('immune_smg.pdf',8,8)
whisker(immune_list, smg_table)
dev.off()
pdf('vascular_smg.pdf',8,8)
whisker(vascular_list, smg_table)
dev.off()
pdf('mes_smg.pdf',8,8)
whisker(mes_list, smg_table)
dev.off()


saveRDS(airway_obj, 'airway_obj.rds')
airway_obj <- readRDS('airway_obj.rds')
#pdf('alveolar.pdf',8,8)
#whisker(alveolar_list, final_table2)
#dev.off()
################################################################################
# Generating Heatmaps for Genes in Airway Epethlial subset across Disease Conditions

basal$Group <- factor(basal$Group, levels = c("HC", "BOS", "CF", 'COPD', 'Asthma'))
markers <- c('ADAM28', 'SOX2', 'EGFR', 'MUC1',
             'TC2N', 'LTB4R','FSCN1','SOX9', 'FOS', 'THBS2', 'ADAMTS1', 'CSPG4', 'SLC2A1', 'CXCL6', 'CP', 'CFB', 'BPIFA1', 'LTF' ,'RARRES1', 'HP',
             'ANO1', 'UPK1B', 'KRT7', 'PCNA', 'SCEL', 'MKI67' )
DoHeatmap(subset(basal, downsample = 30000), features=markers , group.by='Group', assay = 'SCT')+ scale_fill_gradientn(colors = c("darkblue", "lightblue", "red", "darkred"))

#Plot the average expression
basal_average$Group <- factor(basal_average$Group, levels = c("HC", "BOS", "CF", 'COPD', 'Asthma'))
basal_average <- AverageExpression(basal, assays = 'SCT', group.by = 'Group', return.seurat = TRUE)
DoHeatmap(subset(basal_average, downsample = 30000), features=markers , group.by='Group', assay = 'SCT')+ scale_fill_gradientn(colors = c("darkblue", "lightblue", "red", "darkred"))





####-------Heatmaps - Curated-------####

# Create heatmap function
custom_heatmap <- function(cell_type, markers) {
  
  # Subset cell type
  object <- subset(airway_obj, subset = subtype %in% c(cell_type))
  
  # Reorder Levels
  object$Group <- factor(object$Group, levels = c("HC", "BOS", "CF", 'COPD', 'Asthma'))
  
  # Generate first heatmap
  heatmap1 <- DoHeatmap(subset(object, downsample = 30000), markers, group.by='Group', assay = 'SCT', draw.lines = FALSE, raster = FALSE) + 
    theme(axis.text.y = element_text(face="bold", size=12), legend.position = 'right', legend.text = element_text(size=10), 
          legend.title = element_text(size=10), legend.direction = "vertical") + 
    scale_fill_gradientn(colors = c("blue", "white", "red")) + guides(color="none")
  
  print(heatmap1)
  
  # Create Average Expression Object
  object_average <- AverageExpression(object, return.seurat = TRUE, assay = 'SCT', group.by = 'Group')
  
  # Generate second heatmap
  heatmap2 <- DoHeatmap(object_average, markers, draw.lines = FALSE, raster = FALSE, assay = 'SCT') + 
    theme(axis.text.y = element_text(face="bold", size=12), legend.position = 'right', 
          legend.text = element_text(size=10), legend.title = element_text(size=10), 
          legend.direction = "vertical") + scale_fill_gradientn(colors = c("blue", "white", "red")) +
    guides(color="none")
  
  print(heatmap2)
}

# Create DotPlot function
setwd("~/Spatial_Data/Plot1_cell_VlnPlots/Dotplots_All_Groups")
airway_obj$Group <- factor(airway_obj$Group, levels = c("HC", "BOS", "CF", 'COPD', 'Asthma'))
plot1cell_dotplot1 <- function(cell_type, markers) {
  object <- subset(airway_obj, subset = subtype %in% c(cell_type))
  object <- subset(object, subset = Group %ni% c('Asthma'))
  png(filename = paste0(cell_type, '_Dot_Disease.png'), width = 5, height = 10,units = 'in', res = 300)
  plot1cell::complex_dotplot_multiple(object, features = markers, group='Group', color.palette = color.palette)
  dev.off()}
##---Basal---##

# Set markers
markers <- c('ADAM28', 'SOX2', 'EGFR', 'MUC1',
             'TC2N', 'LTB4R','FSCN1','SOX9', 'FOS', 'THBS2', 'ADAMTS1', 'CSPG4', 'SLC2A1', 'CXCL6', 'CP', 'CFB', 'BPIFA1', 'LTF' ,'RARRES1', 'HP',
             'ANO1', 'UPK1B', 'KRT7', 'PCNA', 'SCEL', 'MKI67' )

# Run heatmap function
#custom_heatmap("Basal", markers)

# Run Plot1Cell Violin Plot Function
plot1cell_dotplot1("Basal", markers)

##---Ciliated---###

# Set markers
markers <- c("CYP2F1","MUC1","HMGCS1","CD63","CD9","CCR7","FSCN1","CSPG4","CXCR5",
             "LILRB2","S100B","PLA2G4F","SMIM24","CD163","TM4SF4","MUC5B","CFB"
             ,"BPIFA1","LTF","UPK1B","SCEL","ANO1","KRT7","DMBT1","NTRK2")

# Run heatmap function
#custom_heatmap("Cil", markers)

# Run Plot1Cell Violin Plot Function
plot1cell_dotplot1("Cil", markers)

##---Secretory---###

# Set markers
markers <- c("SERPINA3", "KRT14", "KRT15", "CD8A", 
             "KRT5", "ACTA2", "S100A4","FASN","CD38", "SFTPC", "SFTPD", "BPIFA1","CSTA","HP", 
             "PEBP4", "CXCL5", "MUC5B", "KIT", "EGF", "KDR")

# Run heatmap function
#custom_heatmap("Sec", markers)

# Run Plot1Cell Violin Plot Function
plot1cell_dotplot1("Sec", markers)


##---AM---###

# Set markers
markers <- c("ANPEP", "SYK", "H1F1A", "PIM1", "FOS", "S100A4", "CD14", "FCGR1A", 
             "ATP1B1", "IL33", "LYVE1", "TM4SF18", "N1D1", "GPR183", "AREG", 
             "DMBI1", "CXCL10", "CLEC4F", "VIM", "PLA2G7", "SLC1A3", "FABP3", 
             "CD9", "PDCD1LG2", "CHIT1", "IL1A", "PROX1")

# Run heatmap function
#custom_heatmap("AM", markers)



##---BC---###

# Set markers
markers <- c("TMEM100", "VIM", "S100A12", "SERPINA3", "CD63", "FGF7", "FGFR4", 
             "CCDC78", "FCN3", "KRT7", "MUC1", "CXCL9", "SMIM24", "CD163", 
             "EPCAM", "IGFBP3", "MARCO", "CCR7", "LTF", "TCL1A", "CSTA", 
             "SFTPC", "IL1RL1", "CCN2", "CD1C", "LILRA4", "MYC", "MPEG1", 
             "LYVE1", "TNFRSF17")

# Run heatmap function
#custom_heatmap("BC", markers)

# Run Plot1Cell Violin Plot Function
plot1cell_dotplot1("BC", markers)


##---cDC---###

# Set markers
markers <- c("SERPINA3", "HMGCS1", "CTTN", "FGF7", 
             "SOX2", "KRT5", "FCER1A", "ALOX5AP","CD1A", 
             "KRT14", "TMC5", "CD1C", "MS4A4A", "CXCL6","SOX9","CP", 
             "CTLA4", "PIM2", "GPR171", 
             "AREG", "CXCL10", "GPR183", "MZB1", 
             "PI3", "CCN2", "MYH11", "DES", "SELP", "DMBT1", "GLIPR2", "TGFB3")

# Run heatmap function
#custom_heatmap("cDC", markers)

# Run Plot1Cell Violin Plot Function
plot1cell_dotplot1("cDC", markers)


##---GB---###

# Set markers
markers <- c("SERPINA3", "S100A4", "GLCCI1", "FGF7", "SEC11C", "LAMP2", 
             "TMEM100", "MZB1", "CD9", "THY1", "SMIM24", "FCGR1A", 
             "CXCR5", "CXCL13", "CD79A", "TNFRSF13C", "CD19", "LTF", 
             "DMBT1", "PIGR", "SNAI2", "CYP2F1", "MUC5B", "MYC")


# Run heatmap function
#custom_heatmap("GB", markers)

# Run Plot1Cell Violin Plot Function
plot1cell_dotplot1("GB", markers)

##---Mast---###

# Set markers
markers <- c("RAB27B", "PECAM1", "P2RX1", "SERPINA3", "GNG11", "PIM1", "TMEM100", 
             "S100A12", "HMGCS1", "PTGS1", "FOS", "THY1", "CXCL13", "MZB1", 
             "CD79A", "POU2AF1", "SFTPC", "BPIFA1", "CXCL6", "GPR183", "SNAI2")

# Run heatmap function
#custom_heatmap("Mast", markers)

# Run Plot1Cell Violin Plot Function
plot1cell_dotplot1("Mast", markers)


##---PMN---###

# Set markers
markers <- c("MCEMP1", "ITGAM", "S100A12", "S100A4", "PRKCD", "STEAP4", "MPEG1", 
             "HIF1A", "GPR34", "ELF3", "CD8B", "FOS", "PI3", "STC1", "BPIFA1", 
             "AREG", "ANO1", "FCER1A", "SNAI2", "CD1C", "MMP12")

# Run heatmap function
#custom_heatmap("PMN", markers)

# Run Plot1Cell Violin Plot Function
plot1cell_dotplot1("PMN", markers)



##---Mono/Macs---###

# Set markers
markers <- c("CD14", "VSIG4", "PRKCD", "S100A12", "CD300E", "LTA4H", 
             "HIF1A", "SMIM24", "LAG3", "FOXI1", "CD8B", "MZB1", 
             "CD79A", "POU2AF1", "PI3", "TNFRSF13C", "CCN2", "LILRA4", 
             "GZMB", "CLEC10A", "CD1C", "FCER1A", "CD1A", "GPR171", 
             "MYC", "APOD")

# Run heatmap function
#custom_heatmap("Mono/Macs", markers)

# Run Plot1Cell Violin Plot Function
Mono <- subset(airway_obj, subset = subtype %in% c('Mono/Macs'))
Mono <- subset(Mono, subset = Group %ni% c('Asthma'))
png(filename = paste0('Mono_Macs', '_Dot_Disease.png'), width = 5, height = 10,units = 'in', res = 300)
plot1cell::complex_dotplot_multiple(Mono, features = markers, group='Group', color.palette = color.palette)
dev.off()




##---T Cell---###

# Set markers
markers <- c("S100A4", "VIM", "CD63", "PIM1", "LCK", "SERPINA3", "HIF1A", 
             "LAMP2", "NKG7", "MZB1", "CD79A", "POU2AF1", "CXCL10", "GZMK", 
             "TNFRSF13C", "CCN2", "ANO1", "CD1C", "S100B", "LTF")


# Run heatmap function
#custom_heatmap(c('TkC', 'ThC'), markers)

# Run Plot1Cell Violin Plot Function
Tcell <- subset(airway_obj, subset = subtype %in% c('TkC', 'ThC'))
Tcell <- subset(Tcell, subset = Group %ni% c('Asthma'))
png(filename = paste0('Tcell', '_Dot_Disease.png'), width = 7, height = 10,units = 'in', res = 300)
plot1cell::complex_dotplot_multiple(Tcell, features = markers, group='Group', color.palette = color.palette)
dev.off()



####-------Heatmaps - Top 10-------####

##---Basal---###

# Set markers
markers <- c("ADAM28", "HIF1A", "EHF", "SOX2", "ENAH", "EGFR", "MUC1", 
             "KLF5", "TC2N", "LTB4R", "FSCN1", "SOX9", "FOS", "THBS2", 
             "ADAMTS1", "AIF1", "VIM", "SLC2A1", "KRT14", "CFB", "CXCL6", 
             "CP", "CFB", "PIGR", "BPIFA1", "LTF", "MUC5B", "RARRES1", 
             "SEMA3C", "WASHC1", "ANO1", "CDH1", "UPK1B", "KRT7", "PCNA", 
             "CD9", "KRT5", "CTTN", "HSP90AA1", "GPI", "GPX2")

# Run heatmap function
custom_heatmap("Basal", markers)


##---Ciliated---###

# Set markers
markers <- c("CYP2F1", "FOXJ1", "AGR3", "MUC1", "HMGCS1", "CD63", "TC2N", 
             "EHF", "HSP90AA1", "TSG101", "FOS", "SERPINA3", "JUN", "KRT14", 
             "SLC2A1", "IGFBP3", "MUC5B", "CFB", "CCDC78", "PIGR", "CP", 
             "BPIFA1", "SMPD3", "HSP90AB1", "LTF", "CCN2", "UPK1B", "CD9", 
             "STEAP4", "KLF5", "KRT5", "F3", "CD24", "KRT7", "IL33", 
             "TACSTD2", "F3", "CSTA", "KRT15", "KLK11", "ADAM28", "GPX2")

# Run heatmap function
custom_heatmap("Cil", markers)


##---Secretory---###

# Set markers
markers <- c("PI3", "KRT7", "SERPINA3", "LTF", "KRT14", "TACSTD2", "RARRES1", 
             "CTTN", "KRT5", "VIM", "SFTPC", "MUC1", "SFTPD", "HP", "AGR3", 
             "MALL", "FASN", "PEBP4", "DUOX1", "CXCL6", "CFB", "CP", "CXCL6", 
             "RARRES1", "FOS", "STEAP4", "ANPEP", "JUN", "TMC5", "ELF3", 
             "EHF", "SEC11C", "DMBT1", "LTF", "PIGR", "ANO1", "S100B", 
             "SOX9", "SEC22B", "KLF5", "APOD", "PIGR", "S100B", "BPIFA1")


# Run heatmap function
custom_heatmap("Sec", markers)


##---BC---###

# Set markers
markers <- c("TMEM100", "VIM", "BANK1", "PECAM1", "FCMR", "SERPINA3", "CD63", 
             "MS4A1", "IRF8", "LAMP2", "COL5A2", "THY1", "FBN1", "CCN2", 
             "NID1", "CXCL13", "SFRP2", "IL7R", "AIF1", "LTBP2", "CXCL13", 
             "TNFRSF13C", "HIF1A", "CCR7", "FCMR", "POU2AF1", "CD19", "IRF8", 
             "IL7R", "KLRB1", "SFTPC", "CCN2", "CD1C", "GPR183", "SLIT3", 
             "MYC", "MPEG1", "ENAH", "PECAM1", "CD24", "SPIB", "ALOX5", "APOD")

# Run heatmap function
custom_heatmap("BC", markers)


##---cDC---###

# Set markers
markers <- c("SERPINA3", "VIM", "S100A4", "CD9", "APOD", "SEMA3C", "ITGB4", 
             "HMGCS1", "SEMA3B", "CTTN", "FCER1A", "AREG", "AIF1", "ALOX5AP", 
             "GPR183", "ALOX5", "CD4", "RUNX3", "TACSTD2", "ATP1B1", "S100B", 
             "JUN", "AREG", "SEMA3B", "GPR183", "NID1", "THBS2", "HAVCR2", 
             "IRF8", "ITGB4", "CCN2", "DGKG", "GLIPR2", "S100A4", "SEC22B", 
             "MYO6", "TC2N", "CTNNB1", "APOD", "FGFBP2", "NTRK2")

# Run heatmap function
custom_heatmap("cDC", markers)


##---GB---###

# Set markers
markers <- c("SERPINA3", "APOD", "S100A4", "CD63", "GLCCI1", "LTF", "PLA2G2A", 
             "FGF7", "DMBT1", "PIGR", "COL5A2", "THY1", "FBN1", "CD163", 
             "NID1", "SFRP2", "CD79A", "POU2AF1", "DNAJB9", "TNFRSF13C", 
             "PECAM1", "VIM", "AREG", "MZB1", "SFTPC", "CD27", "DMBT1", 
             "LTF", "TNFRSF17", "PIGR", "SERPINA3", "CD63", "SVEP1", "IGF1", 
             "ADAM28", "MYC", "APOD", "FCMR", "COL5A2")

# Run heatmap function
custom_heatmap("GB", markers)


##---Mono/Macs---###

# Set markers
markers <- c("RAB27B", "PECAM1", "VIM", "P2RX1", "S100A4", "SERPINA3", "GNG11", 
             "PIM1", "CTTN", "TMEM100", "CCN2", "FBN1", "MS4A4A", "PDCD1LG2", 
             "PIM2", "AREG", "KIT", "IL1RL1", "RAB44", "SFTPC", "WASHC1", 
             "LTC4S", "SEC11C", "P2RX1", "ALOX5AP", "MS4A2", "VIM", "FCER1A", 
             "GPR183", "CCN2", "SLC18A2", "HPGDS", "KIT", "S100A4", "CD274")

# Run heatmap function
custom_heatmap("Mono/Macs", markers)


##---Mast---###

# Set markers
markers <- c("S100A4", "VIM", "PECAM1", "CD14", "FCGR3A", "CLEC12A", "VSIG4", 
             "ITGAM", "PRKCD", "CD63", "AREG", "CCN2", "GPR183", "MYC", 
             "CLEC4E", "AREG", "GPR183", "SLC1A3", "CLEC10A", "ADAM28", 
             "HIF1A", "FOS", "CD4", "LYVE1", "CCN2", "CD1C", "FCER1A", 
             "CLEC10A", "PTGS1", "GPR183", "SLIT3", "SLC1A3", "CD4", "MYC")

# Run heatmap function
custom_heatmap("Mast", markers)


##---PMN---###

# Set markers
markers <- c("MCEMP1", "VIM", "ALOX5", "CLEC12A", "ANPEP", "ITGAM", "S100A12", 
             "S100A4", "PRKCD", "TMEM100", "CLEC4E", "CCN2", "CD163", "FOS", 
             "MPEG1", "APOLD1", "IL33", "FCGR3A", "IL7R", "LILRB2", "CLEC4E", 
             "APOLD1", "FOS", "AQP9", "HIF1A", "MTUS1", "RALB", "IL33", 
             "SDCBP", "VWF", "CCN2", "SLIT3", "CD4", "MYC", "LTBP2", "IL33", 
             "SVEP1", "CLEC4E", "APOLD1", "SFTPC")

# Run heatmap function
custom_heatmap("PMN", markers)


##---T Cell---###

# Set markers
markers <- c("S100A4", "VIM", "PECAM1", "CD63", "TMEM100", "PIM1", "LCK", 
             "SERPINA3", "FCN3", "SNX3", "KLRD1", "GZMB", "FCGR3A", "CCN2", 
             "AREG", "FBN1", "GPR171", "CTLA4", "FCMR", "GPR183", "IL7R", 
             "SFRP2", "MZB1", "GLCCI1", "CXCR6", "GZMK", "CCN2", "SFTPC", 
             "AGER", "IL7R", "FGFBP2", "SLIT3", "SVEP1", "IL33", "TC2N", 
             "GZMB", "APOD")

# Run heatmap function
custom_heatmap(c('TkC', 'ThC'), markers)


####-------Control vs CF-------####

# Create New Ident column and set to Idents #
airway_obj$celltype.sample <- paste(airway_obj$subtype, airway_obj$Group, sep = "_")
Idents(airway_obj) <- "celltype.sample"

# Rename specific cell types
airway_obj$celltype.sample <- gsub("TkC_HC", "T_HC", airway_obj$celltype.sample)
airway_obj$celltype.sample <- gsub("ThC_HC", "T_HC", airway_obj$celltype.sample)
airway_obj$celltype.sample <- gsub("TkC_CF", "T_CF", airway_obj$celltype.sample)
airway_obj$celltype.sample <- gsub("ThC_CF", "T_CF", airway_obj$celltype.sample)

# Make Cell and Stim List #
cf_list <- c("Basal_CF", "Cil_CF", "Sec_CF", "BC_CF", "cDC_CF", "Mono/Macs_CF", 
             "Mast_CF", "PMN_CF", "T_CF", "GB_CF")
ctrl_list <- c("Basal_HC", "Cil_HC", "Sec_HC", "BC_HC", "cDC_HC", "Mono/Macs_HC", 
               "Mast_HC", "PMN_HC", "T_HC", "GB_HC") 

# Subset airway object for only desired groups
CF_O <- subset(airway_obj, subset = Group %in% c("CF", "HC"))

# PrepSCTFindMarkers # 
CF_O <- PrepSCTFindMarkers(CF_O)

# Set idents
Idents(CF_O) <- "celltype.sample"

# Set Assay #
DefaultAssay(CF_O) <- "SCT"

# Run FindMarkers CF v Control 
for (i in 1:10) {
  name = FindMarkers(CF_O, ident.1 = cf_list[[i]], ident.2 = ctrl_list[[i]], verbose = TRUE)
  file_name = gsub("/", "_", unlist(strsplit(cf_list[[i]], split = "_"))[1])
  file = paste0(file_name, "_wilcox_", "_CF_", ".csv")
  write.csv(name, file)
}

# Run FindMarkers Control v CF 
for (i in 1:10) {
  name = FindMarkers(CF_O, ident.1 = ctrl_list[[i]], ident.2 = cf_list[[i]], verbose = TRUE)
  file_name = gsub("/", "_", unlist(strsplit(cf_list[[i]], split = "_"))[1])
  file = paste0(file_name, "_wilcox_", "_Ctrl_", ".csv")
  write.csv(name, file)
}

# Create heatmap function
custom_heatmapv2 <- function(cell_type, markers) {
  
  # Subset cell type
  object <- subset(CF_O, subset = subtype %in% c(cell_type))
  
  # Reorder Levels
  levels(object) <- c("HC", "CF")
  
  # Generate first heatmap
  heatmap1 <- DoHeatmap(subset(object, downsample = 15000), markers, assay = 'SCT', draw.lines = FALSE, raster = FALSE) + 
    theme(axis.text.y = element_text(face="bold", size=12), legend.position = 'right', legend.text = element_text(size=10), 
          legend.title = element_text(size=10), legend.direction = "vertical") + 
    scale_fill_gradientn(colors = c("blue", "white", "red")) + guides(color="none")
  
  print(heatmap1)
  
  # Create Average Expression Object
  object_average <- AverageExpression(object, return.seurat = TRUE, assay = 'SCT')
  
  # Generate second heatmap
  heatmap2 <- DoHeatmap(object_average, markers, draw.lines = FALSE, raster = FALSE, assay = 'SCT') + 
    theme(axis.text.y = element_text(face="bold", size=12), legend.position = 'right', 
          legend.text = element_text(size=10), legend.title = element_text(size=10), 
          legend.direction = "vertical") + scale_fill_gradientn(colors = c("blue", "white", "red")) +
    guides(color="none")
  
  print(heatmap2)
}

# Create markers function
marker_v1 <- function(file_path) {
  # Read the data from the file
  data <- read_csv(file_path)
  top_genes <- data %>%
    filter(pct.1 >= 0.1, avg_log2FC > 0) %>%
    arrange(p_val_adj) %>%
    head(10) %>%
    pull(...1)  
  return(as.character(top_genes))
}

####-------Heatmaps - CF vs Control-------####

##---Basal---###

# Set markers
markers <- c(marker_v1("DE_CF/Ctrl/Basal_wilcox__Ctrl_.csv"), 
             marker_v1("DE_CF/CF/Basal_wilcox__CF_.csv"))

# Run heatmap function
custom_heatmapv2("Basal", markers)


##---Ciliated---###

# Set markers
markers <- c(marker_v1("DE_CF/Ctrl/Cil_wilcox__Ctrl_.csv"), 
             marker_v1("DE_CF/CF/Cil_wilcox__CF_.csv"))

# Run heatmap function
custom_heatmapv2("Cil", markers)


##---Secretory---###

markers <- c(marker_v1("DE_CF/Ctrl/Sec_wilcox__Ctrl_.csv"), 
             marker_v1("DE_CF/CF/Sec_wilcox__CF_.csv"))

# Run heatmap function
custom_heatmapv2("Sec", markers)


##---BC---###

markers <- c(marker_v1("DE_CF/Ctrl/BC_wilcox__Ctrl_.csv"), 
             marker_v1("DE_CF/CF/BC_wilcox__CF_.csv"))

# Run heatmap function
custom_heatmapv2("BC", markers)


##---cDC---###

# Set markers
markers <- c(marker_v1("DE_CF/Ctrl/cDC_wilcox__Ctrl_.csv"), 
             marker_v1("DE_CF/CF/cDC_wilcox__CF_.csv"))

# Run heatmap function
custom_heatmapv2("cDC", markers)


##---GB---###

# Set markers
markers <- c(marker_v1("DE_CF/Ctrl/GB_wilcox__Ctrl_.csv"), 
             marker_v1("DE_CF/CF/GB_wilcox__CF_.csv"))

# Run heatmap function
custom_heatmapv2("GB", markers)

##---Mono/Macs---###

# Set markers
markers <- c(marker_v1("DE_CF/Ctrl/Mono_Macs_wilcox__Ctrl_.csv"), 
             marker_v1("DE_CF/CF/Mono_Macs_wilcox__CF_.csv"))

# Run heatmap function
custom_heatmapv2("Mono/Macs", markers)


##---Mast---###

# Set markers
markers <- c(marker_v1("DE_CF/Ctrl/Mast_wilcox__Ctrl_.csv"), 
             marker_v1("DE_CF/CF/Mast_wilcox__CF_.csv"))

# Run heatmap function
custom_heatmapv2("Mast", markers)


##---PMN---###

# Set markers
markers <- c(marker_v1("DE_CF/Ctrl/PMN_wilcox__Ctrl_.csv"), 
             marker_v1("DE_CF/CF/PMN_wilcox__CF_.csv"))

# Run heatmap function
custom_heatmapv2("PMN", markers)


##---T Cell---###

# Set markers
markers <- c(marker_v1("DE_CF/Ctrl/T_wilcox__Ctrl_.csv"), 
             marker_v1("DE_CF/CF/T_wilcox__CF_.csv"))

# Run heatmap function
custom_heatmapv2(c('TkC', 'ThC'), markers)


####-------Asthma Issue-------####

# Create individual objects based on disease type
HC <- subset(airway_obj, subset = Group %in% "HC")
BOS <- subset(airway_obj, subset = Group %in% "BOS")
CF <- subset(airway_obj, subset = Group %in% "CF")
COPD <- subset(airway_obj, subset = Group %in% "COPD")
Asthma <- subset(airway_obj, subset = Group %in% "Asthma")

# Create List of Seurat Objects #
Disease_List <- c(HC, BOS, CF, COPD, Asthma)

# Total sequencing depth function #
my_function <- function(x) {
  print(unique(x$ID))
  print(x)
  print(sum(x$nCount_SCT))
  print(mean(x$nCount_SCT))
}

# Apply to list #
lapply(Disease_List, my_function)


####-------Dotplots - CF vs Control-------####
setwd("~/Spatial_Data/Plot1_cell_VlnPlots/Dotplot_CF_HC")
plot1cell_dotplot2 <- function(cell_type, markers) {
  object <- subset(airway_obj, subset = subtype %in% c(cell_type))
  object <- subset(object, subset = Group %ni% c('Asthma'))
  object$Group <- factor(object$Group, levels = c("HC", "CF", 'COPD', 'BOS'))
  png(filename = paste0(cell_type, '_Dot_CF_HC.png'), width = 5, height = 10,units = 'in', res = 300)
  plot1cell::complex_dotplot_multiple(object, features = markers, group='Group', color.palette = color.palette)
  dev.off()}
##---Basal---###

# Set markers
markers <- c('CD34', 'ADAM28', 'FBN1', 'LGR6', 'GPX2', 'FOS', 'LTBP2', 'CYSLTR1', 'DMBT1', 'SOX9', 'VSIG4', 'CYP2F1', 'FOXJ1', 'KRT5', 'AGR3' ,
             'S100A7', 'HP', 'CXCL5', 'MS4A2', 'CXCL10', 'IL1RL1', 'KIT' ,'AREG', 'LTF', 'HPGDS', 'STC1', 'CXCL9', 'NTRK2', 'BPIFA1', 'SFRP2')


# Run Dotplot function
plot1cell_dotplot2("Basal", markers)


##---Cil---###

# Set markers
markers <- c('ANPEP', 'KRT5', 'FOS', 'CYP2F1', 'ADAM28', 'KRT15', 'PI3', 'HMGCS1', 'GPX2', 'DMBT1', 'S100A4', 'SLC7A11', 'POU2AF1', 'PRKCA', 'IL33A1-exon1A1-exon2',
             'MZB1', 'LTF', 'VIM', 'BPIFA1', 'ACE', 'SPDEF', 'CXCL5', 'STAT4', 'RAB27B', 'SEMA3B', 'GSDMB', 'TFAP2C', 'CFB', 'MUC5B', 'AIF1' )


# Run Dotplot function
plot1cell_dotplot2("Cil", markers)

##---Sec---###

# Set markers
markers <- c('MYH11' ,'LMOD1', 'CNN1', 'LGR6', 'ACTA2', 'DMBT1', 'APOD', 'S100B', 'LYVE1', 'P2RX1', 'NTRK2', 'SERPINA3', 'FGF7', 'ANO1', 'KLK11',
             'CXCL5', 'SFTPC', 'S100A7', 'FCER1A', 'RND1', 'SFTPD', 'SPIB', 'MCEMP1', 'SELL', 'NCEH1', 'CLEC4E', 'RGS5', 'AQP9', 'CD79A', 'CXCL6' )


# Run Dotplot function
plot1cell_dotplot2("Sec", markers)

##---Mast---###

# Set markers
markers <- c('ASCL3', 'FOXI1', 'SERPINA3', 'PLA2G2A', 'VSIG4', 'MARCO', 'S100A12', 'PDGFRB', 'RAMP2', 'FGF7', 'PDGFRA', 'CFTR', 'ECSCR', 'RAB27B', 'PIM1',
             'CXCL10', 'CD79A', 'SPDEF', 'UPK1B', 'SFTPC', 'HP', 'BANK1', 'PI3', 'BPIFA1', 'GZMK', 'POU2AF1', 'IL7R', 'CD86', 'MKI67', 'GSDME') 

# Run Dotplot function
plot1cell_dotplot2("Mast", markers)

##---BC---###

# Set markers
markers <- c('GPX2', 'IL13', 'FOXJ1', 'KRT5', 'ACKR1', 'FGFR4', 'TM4SF18', 'PLA2G2A', 'GZMB', 'SERPINA3', 'TMEM100', 'NTRK2', 'MCEMP1', 'KDR', 'VSIG4',
             'CXCL13', 'TCL1A', 'MMP9', 'LTF', 'IL2RA', 'P2RX1', 'SFTPC', 'CD1A', 'ARL14', 'CD38', 'KIT', 'CSTA', 'CCR7', 'MKI67', 'GPR34' ) 

# Run Dotplot function
plot1cell_dotplot2("BC", markers)

##---GB---###

# Set markers
markers <- c('SERPINA3', 'TREM2', 'VSIG4', 'KRT14', 'KRT5', 'DMBT1', 'PCOLCE2', 'RAMP2', 'ICA1', 'ACKR1', 'NTRK2', 'KRT15', 'STEAP4', 'TMEM100', 'MMRN1',
             'CCR7', 'CP', 'BPIFA1', 'CLEC4E', 'CXCL5', 'SFTPC', 'IGFBP3', 'CD247', 'TNFRSF13C', 'CXCL6', 'BANK1', 'CTLA4', 'THY1', 'IL7R', 'CXCL10'  ) 

# Run Dotplot function
plot1cell_dotplot2("GB", markers)

##---cDC---###

# Set markers
markers <- c('LGR5', 'SERPINA3', 'TMEM100', 'SOX2', 'CNN1', 'CDH1', 'PRKCA', 'LTF', 'CTSL', 'ITGB4', 'S100A4', 'MET','FGF7', 'NTRK2', 'SEMA3C',
             'AREG', 'FCER1A', 'POU2AF1', 'CTLA4', 'CD1C', 'CXCL10', 'ANO1', 'KIT', 'CCR7', 'MZB1', 'ACE', 'ELF3', 'CD2', 'GSDMC', 'TNFRSF13C' ) 

# Run Dotplot function
plot1cell_dotplot2("cDC", markers)

##---PMN---###

# Set markers
markers <- c('GPX2', 'FCGR1A', 'CD300E', 'SELE', 'CD14', 'SERPINA3', 'STEAP4', 'VSIG4', 'PIM1', 'S100A4', 'CYP2F1', 'SEMA3B', 'KRT15', 'RARRES1', 'AGR3',
             'S100A7', 'CXCL13', 'CD1A', 'GJA5', 'STC1', 'APOLD1', 'MMP12', 'SLC7A11', 'GKN2', 'CXCL5', 'BPIFA1', 'CD79A', 'CD19', 'PROX1', 'LILRA4' ) 

# Run Dotplot function
plot1cell_dotplot2("PMN", markers)

##---Mono/Macs---###

# Set markers
markers <- c('SERPINA3', 'KRT5', 'PIM1', 'CTSL', 'VSIG4', 'STEAP4', 'KRT15', 'FOXJ1', 'CD14', 'TREM2', 'KLF5', 'AGR3', 'FCGR1A', 'ETV5', 'MARCO',
             'IL2RA', 'CD79A', 'POU2AF1', 'CD27', 'AREG', 'MZB1', 'CD19', 'CLEC4E', 'HP', 'CTLA4', 'CCR7' ,'TNFRSF17', 'PI3', 'GPR171', 'CXCL10' ) 

# Run Dotplot function
Mono <- subset(airway_obj, subset = subtype %in% c('Mono/Macs'))
Mono <- subset(Mono, subset = Group %ni% c('Asthma'))
Mono$Group <- factor(Mono$Group, levels = c("HC", "CF", 'COPD', 'BOS'))
png(filename = paste0('Mono_Macs', '_Dot_CF_HC.png'), width = 5, height = 10,units = 'in', res = 300)
plot1cell::complex_dotplot_multiple(Mono, features = markers, group='Group', color.palette = color.palette)
dev.off()

##---Tcell---###

# Set markers
markers <- c('SERPINA3', 'GPX2', 'KRT5', 'KLK11', 'KRT15', 'VSIG4', 'PLA2G2A', 'KLRC1', 'S100A12', 'CYP2F1', 'FASLG', 'PCOLCE2', 'KLF5', 'IL33FULL-exon7-exon8', 'TP63',
             'CXCL13', 'CD79A', 'CXCL10', 'POU2AF1', 'SFTPC', 'SPIB', 'MMP9', 'MZB1', 'SELL', 'TNFRSF13C', 'CD19', 'CLEC4E', 'TNFRSF13B', 'CCR7', 'TNFRSF17' ) 

# Run Dotplot function
Tcell <- subset(airway_obj, subset = subtype %in% c('ThC', 'TkC'))
Tcell <- subset(Tcell, subset = Group %ni% c('Asthma'))
Tcell$Group <- factor(Tcell$Group, levels = c("HC", "CF", 'COPD', 'BOS'))
png(filename = paste0('TCell', '_Dot_CF_HC.png'), width = 7, height = 10,units = 'in', res = 300)
plot1cell::complex_dotplot_multiple(Tcell, features = markers, group='Group', color.palette = color.palette)
dev.off()

####### Single Cell RL Validation #############################################
copd_hc <- readRDS("/media/haider/Files/Spatial_Data/GSE173896_COPD.rds")

# Separate Basal Cells
FeaturePlot(copd_hc,c('TP63'),  cols = c('grey94', 'darkblue'), order=T)
plot <- DimPlot(copd_hc, label=TRUE)
copd_hc <- CellSelector(plot = plot, object = copd_hc, ident = 'Basal')

# Subset for Airway cells
copd_hc$celltypes <- copd_hc@active.ident
copd_airway <- subset(copd_hc, celltypes %in% c('Basal', 'Macrophage-1','Ciliated' ,'Macrophage-2', 'cMonocyte','Basal/Club', 'MAST', 'Bcell','Basophil',   "Nerutrophil" ))
saveRDS(copd_airway, "~/Squidpy_Spatial_analysis/Airway/Cellchat/copd_airway.rds")
#Basal, Basal/Club, Cil, MAST, MAcrhophage 1, Mac 2, cmonocyte, Bcell, Basophil, Neutrophil


# Create Individual CellChat objects for Ctrl
ctrl_lis <- subset(copd_airway, Class %in% c("Non-smoker"))
data.inpu1 <- GetAssayData(ctrl_lis, assay='RNA', slot='data')
labels1=Idents(ctrl_lis)
meta1 <- data.frame(group=labels1, row.names=names(labels1))
ctrl_chat <- createCellChat(object=data.inpu1, meta=meta1, group.by = 'group')
# CellChatDB for cell-cell communication analysis
# use all CellChatDB for cell-cell communication analysis
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

# set the used database in the object
ctrl_chat@DB <- CellChatDB.use

ctrl_chat <- subsetData(ctrl_chat) # subset the expression data of signaling genes for saving computation cost
plan("multisession", workers=10)

ctrl_chat <- identifyOverExpressedGenes(ctrl_chat)
ctrl_chat <- identifyOverExpressedInteractions(ctrl_chat)
# Project gene expression data onto protein-protein interactions (PPI-human)
ctrl_chat <- smoothData(ctrl_chat, adj=PPI.human)

# Compute the communication probability and infer cell communication network
ctrl_chat <- computeCommunProb(ctrl_chat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#ctrl_chat <- filterCommunication(ctrl_chat, min.cells = 10)
#Infer cell-cell communication at signalling pathway level
ctrl_chat <- computeCommunProbPathway(ctrl_chat)
# Calculate the aggregate cell-cell communication network
ctrl_chat <- aggregateNet(ctrl_chat)
ctrl_chat

################################################################################
# Create Individual CellChat objects for copd
copd_lis <- subset(copd_airway, subset=Class %in% 'COPD')
data.inpu2 <- GetAssayData(copd_lis, assay='RNA', slot='data')
labels2=Idents(copd_lis)
meta2 <- data.frame(group=labels2, row.names=names(labels2))
copd_chat <- createCellChat(object=data.inpu2, meta=meta2, group.by = 'group')
# CellChatDB for cell-cell communication analysis
# use all CellChatDB for cell-cell communication analysis
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB

# set the used database in the object
copd_chat@DB <- CellChatDB.use

copd_chat <- subsetData(copd_chat) # subset the expression data of signaling genes for saving computation cost
plan("multisession", workers=10)

copd_chat <- identifyOverExpressedGenes(copd_chat)
copd_chat <- identifyOverExpressedInteractions(copd_chat)
# Project gene expression data onto protein-protein interactions (PPI-mouse)
copd_chat <- smoothData(copd_chat, adj=PPI.human)

# Compute the communication probability and infer cell communication network
copd_chat <- computeCommunProb(copd_chat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
#copd_chat <- filterCommunication(copd_chat, min.cells = 10)
#Infer cell-cell communication at signalling pathway level
copd_chat <- computeCommunProbPathway(copd_chat)
# Calculate the aggregate cell-cell communication network
copd_chat <- aggregateNet(copd_chat)
################################################################################
# Merge the both Cellchat objects into one single Cellchat object
object.list <- list(ctrl=ctrl_chat, copd=copd_chat)

cellchat
saveRDS(cellchat, 'copd_hc_cellchat_subset.rds')
saveRDS(object.list, file = "cellchat_object.list_subset.rds")
cellchat <- readRDS("~/Squidpy_Spatial_analysis/Airway/Cellchat/copd_hc_cellchat_subset.rds")
object.list <- readRDS("~/Squidpy_Spatial_analysis/Airway/Cellchat/cellchat_object.list_subset.rds")

group.celltype <- c("Basal","Nerutrophil", "Macrophage-1", "MAST","Macrophage-2", "cMonocyte","Ciliated","Basal/Club",
                    "Bcell","Basophil") 
group.celltype <- factor(group.celltype, levels = c("Basal","Nerutrophil", "Macrophage-1", "MAST","Macrophage-2", "cMonocyte","Ciliated","Basal/Club",
                                                    "Bcell","Basophil"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.celltype)})
cellchat <- mergeCellChat(object.list, add.names=names(object.list))
# Circle Plots
weight.max <- getMaxWeight(object.list, slot.name = c("idents", "net", "net"), attribute = c("idents","count", "count.merged"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count.merged, weight.scale = T, label.edge= T, edge.weight.max = weight.max[3], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

##
ctrl_chat <- object.list[[1]]
copd_chat <-object.list[[2]]
ctrl_chat <-subsetCellChat(ctrl_chat, idents.use = c('Basal', 'MAST'))
copd_chat <-subsetCellChat(copd_chat, idents.use = c('Basal', 'MAST'))
object.list <- list(ctrl=ctrl_chat, copd=copd_chat)
cellchat_bm <- mergeCellChat(object.list, add.names=names(object.list))

# Chord diagram
ctrl_pathways.show <-  cellchat@netP$ctrl$pathways
copd_pathways.show <-  cellchat@netP$copd$pathways
pathways.show <-  copd_pathways.show[copd_pathways.show %in% ctrl_pathways.show]


chord_plot <- function(list, pathway) {
  par(mfrow = c(1,1), xpd=TRUE)
  for (i in 1:length(pathway)) {
    pairLR1 <- extractEnrichedLR(list[[1]], signaling = pathway[i], geneLR.return = FALSE)
    pairLR1 <- pairLR1$interaction_name
    pairLR2 <- extractEnrichedLR(list[[2]], signaling = pathway[i], geneLR.return = FALSE)
    pairLR2 <- pairLR2$interaction_name
    pathway.use =pathway[i]
    for(j in 1:length(pairLR)){
      print(pairLR[j])
      #pdf(file=paste0(pairLR[j], "_chord_diagram.pdf"), width=8, height=8)
      netVisual_individual(obj, signaling =pathway.use, pairLR.use = pairLR[j], layout = "chord")
      #dev.off()
    }}}

par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Basal to Mast
pdf(file=paste0('basal_to_mast', "_chord_diagram.pdf"), width=10, height=10)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 1, targets.use = c(2), lab.cex = 0.5, title.name = paste0("Signaling from Basal - ", names(object.list)[i]))}
dev.off()
par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Mast to Basal
pdf(file=paste0('mast_to_basal', "_chord_diagram.pdf"), width=10, height=10)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 2, targets.use = c(1), lab.cex = 0.5, title.name = paste0("Signaling from Mast - ", names(object.list)[i]))}
dev.off()


setwd("~/Squidpy_Spatial_analysis/Airway/Cellchat/COPD")
chord_plot(copd_chat, pathways.show)
setwd("~/Squidpy_Spatial_analysis/Airway/Cellchat/HC")
chord_plot(ctrl_chat, ctrl_pathways.show)
setwd("/home/Nico/C0_Cord_Plots")
chord_plot(C0_chat, C0_pathways.show)
setwd("/home/Nico/R120_Cord_Plots")
chord_plot(R120_chat, R120_pathways.show)


#Heatmaps


for (j in 1:length(pathways.show)){
  weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show[j]) # control the edge weights across different datasets
  par(mfrow = c(1,2), xpd=TRUE)
  for (i in 1:length(object.list)) {
    netVisual_aggregate(object.list[[i]], signaling = pathways.show[j], layout = "chord", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste(pathways.show[j], names(object.list)[i]))}}
    
    #ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show[j], color.heatmap = "Reds",title.name = paste(pathways.show[j], "signaling ",names(object.list)[i]))}
    #setwd("~/Squidpy_Spatial_analysis/Airway/Cellchat/Heatmaps")
    #pdf(file=paste0('copd_ctrl_',pathways.show[j] , "_heatmap.pdf"), width=15, height=8)
    #ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
    #dev.off()
    