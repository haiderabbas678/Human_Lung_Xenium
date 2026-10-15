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
use_python('/home/haider/.venv/bin/python')
################################################################################p
# Load complete spatial object
final_obj <- readRDS("/media/haider/Files/Spatial_Data/annotated_comb_final.rds")
# Remove low quality clusters 
final_obj <- subset(final_obj, subset=subtype %ni% c('LQ',  "AT2_EC_DB" ))
final_obj <- subset(final_obj, subset= ID %ni% c('MTS43D', 'MTS38D'))
FB_sub <- subset(final_obj, subset = subtype %in% c('FB', 'APOD+_FB', 'PLA2G2A+_FB'))

DimPlot(final_obj, label=TRUE)
FeaturePlot(final_obj, cols = c('grey94', 'darkblue'), features = 'LTBP2')

DotPlot(final_obj, features = 'LTBP2', cols = c('red', 'blue', 'yellow', 'black', 'magenta'))

VlnPlot(FB_sub, features = 'LTBP2', split.by = 'Group', raster = FALSE, pt.size = 0)

levels(factor(final_obj$ID))


files <- c("BOS11A", "BOS12D", "BOS13E" ,"C48A" , "C50C","C51E","C52E","C56E", "ID6","MTS12D","MTS21A","MTS27E","MTS29E","MTS48B","MTS7D" )

# DimPlot for Celltypes
for (x in 1:length(files)) {
  pdf(paste0(files[x],'_spatial_FB_image'),5,5)
  FB_sub@active.ident <- FB_sub$subtype
  sample <- subset(FB_sub, subset=ID %in% files[x])
  print(ImageDimPlot(sample, cols = c('red', 'blue', 'green'), size = 1))
  dev.off()
} 

# FeaturePlot for LTBP2
for (x in 1:length(files)) {
  pdf(paste0(files[x],'_spatial_FB_image_LTBP2'),5,5)
  sample <- subset(FB_sub, subset=ID %in% files[x])
  print(ImageFeaturePlot(sample,features = 'LTBP2', size = 2, cols = c('grey94', 'red')))
  dev.off()
} 
conda_list()
conda_list()[[1]][2] %>% 
  use_condaenv(required = TRUE)
use_condaenv('base')
# Add layers
data <- read_h5ad("~/Squidpy_Spatial_analysis/Airway_subclass_obj")
ep <- data$obs[data$obs$Layers %in% 'Airway_Ep',]
sub <- data$obs[data$obs$Layers %in% 'Sub_Ep',]
smg <- data$obs[data$obs$Layers %in% 'SMG',]

FB_sub$Layers <- sapply(strsplit(as.character(colnames(FB_sub)), "_"), `[`, 2)
FB_sub$Layers[FB_sub$Layers %in% ep$cell_ids] <- 'Airway_Ep'
FB_sub$Layers[FB_sub$Layers %in% sub$cell_ids] <- 'Sub_Ep'
FB_sub$Layers[FB_sub$Layers %in% smg$cell_ids] <- 'SMG'
FB_sub$Layers[FB_sub$Layers %ni% c('Airway_Ep', 'Sub_Ep', 'SMG')] <- 'Non-Airway'

final_obj$Layers <- sapply(strsplit(as.character(colnames(final_obj)), "_"), `[`, 2)
final_obj$Layers[final_obj$Layers %in% ep$cell_ids] <- 'Airway_Ep'
final_obj$Layers[final_obj$Layers %in% sub$cell_ids] <- 'Sub_Ep'
final_obj$Layers[final_obj$Layers %in% smg$cell_ids] <- 'SMG'
final_obj$Layers[final_obj$Layers %ni% c('Airway_Ep', 'Sub_Ep', 'SMG')] <- 'Non-Airway'

# Dimplot for Layers
for (x in 1:length(files)) {
  pdf(paste0(files[x],'_spatial_FB_image_Layers'),5,5)
  FB_sub@active.ident <- factor(FB_sub$Layers)
  sample <- subset(FB_sub, subset=ID %in% files[x])
  print(ImageDimPlot(sample, cols = c('red', 'grey', 'blue', 'green'), size = 1))
  dev.off()
} 

reticulate::install_miniconda()
anndata::install_anndata()
# VlnPlots per Sample for LTBP2 expression
for (x in 1:length(files)) {
  pdf(paste0(files[x],'_spatial_FB_LTBP2_vlnplt'),5,5)
  FB_sub@active.ident <- factor(FB_sub$subtype)
  sample <- subset(FB_sub, subset=ID %in% files[x])
  print(VlnPlot(sample, pt.size = 0, features = 'LTBP2'))
  dev.off()
} 



FB_sub@active.ident <- factor(FB_sub$subtype)
VlnPlot(FB_sub, features = 'LTBP2', pt.size = 0)
create_table <- function(final_table2, total_table){
  final_table2[,1] <- final_table2[,1]/sum(total_table[,1])*100
  final_table2[,2] <- final_table2[,2]/sum(total_table[,2])*100
  final_table2[,3] <- final_table2[,3]/sum(total_table[,3])*100
  final_table2[,4] <- final_table2[,4]/sum(total_table[,4])*100
  final_table2[,5] <- final_table2[,5]/sum(total_table[,5])*100
  final_table2[,6] <- final_table2[,6]/sum(total_table[,6])*100
  final_table2[,7] <- final_table2[,7]/sum(total_table[,7])*100
  final_table2[,8] <- final_table2[,8]/sum(total_table[,8])*100
  final_table2[,9] <- final_table2[,9]/sum(total_table[,9])*100
  final_table2[,10] <- final_table2[,10]/sum(total_table[,10])*100
  final_table2[,11] <- final_table2[,11]/sum(total_table[,11])*100
  final_table2[,12] <- final_table2[,12]/sum(total_table[,12])*100
  final_table2[,13] <- final_table2[,13]/sum(total_table[,13])*100
  final_table2[,14] <- final_table2[,14]/sum(total_table[,14])*100
  final_table2[,15] <- final_table2[,15]/sum(total_table[,15])*100
  final_table2[,16] <- final_table2[,16]/sum(total_table[,16])*100
  final_table2[,17] <- final_table2[,17]/sum(total_table[,17])*100
  final_table2[,18] <- final_table2[,18]/sum(total_table[,18])*100
  final_table2[,19] <- final_table2[,19]/sum(total_table[,19])*100
  final_table2[,20] <- final_table2[,20]/sum(total_table[,20])*100
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
FB_table <- create_table(table(FB_sub$subtype, FB_sub$ID), table(final_obj$subtype, final_obj$ID))


whisker <- function(sample_list, obj){
  for (i in 1:length(sample_list)){
    print(ggplot(obj[obj$Var1 %in% c(sample_list[i]),], aes(x=`Var2`, y=Freq)) +
            geom_boxplot(outlier.shape = NA)+ geom_jitter(width = 0.2)+ggtitle(sample_list[i]))
    }}

pdf('FB_Box_plots_total_cells', 5,5)
whisker(c('APOD+_FB', 'FB', 'PLA2G2A+_FB'), FB_table)
dev.off()

# Epithelial + Sub-Epithelial
final_ep_sub <- subset(final_obj, subset=Layers %in% c('Airway_Ep', 'Sub_Ep'))
FB_ep_sub <- subset(FB_sub, subset=Layers %in% c('Airway_Ep', 'Sub_Ep'))

FB_sub_ep_table <- create_table(table(FB_ep_sub$subtype, FB_ep_sub$ID), table(final_ep_sub$subtype, final_ep_sub$ID))
pdf('FB_Box_plots_total_cells_ep_sub', 5,5)
whisker(c('APOD+_FB', 'FB', 'PLA2G2A+_FB'), FB_sub_ep_table)
dev.off()

# VlnPlots per Sample for LTBP2 expression
for (x in 1:length(files)) {
  pdf(paste0(files[x],'_spatial_FB_LTBP2_sub_ep_vlnplt'),5,5)
  FB_ep_sub@active.ident <- factor(FB_ep_sub$subtype)
  sample <- subset(FB_ep_sub, subset=ID %in% files[x])
  print(VlnPlot(sample, pt.size = 0, features = 'LTBP2'))
  dev.off()}
  
# 


FB_ep_sub@meta.data$Group <- factor(FB_ep_sub@meta.data$Group,  levels = c('HC', 'Asthma', 'BOS', 'CF', 'COPD'))
FB_ep_sub@active.ident <- factor(FB_ep_sub$Group)

VlnPlot(FB_ep_sub, pt.size = 0, features = 'nCount_Xenium')

mean_expression <- AverageExpression(FB_ep_sub, 
                                     # features = "LTBP2", 
                                     group.by = "nCount_SCT")$SCT

# Convert the mean_expression table to a data frame suitable for plotting
mean_df <- data.frame(
  group = colnames(mean_expression),
  mean_val = as.numeric(mean_expression)
)


head(FB_ep_sub$nCount_SCT)


# Generate the violin plot using VlnPlot
  
# nCount_SCT
FB_ep_sub@active.ident <- factor(FB_ep_sub$ID)
p <- VlnPlot2(FB_ep_sub, features = 'nCount_SCT', pt.size = 0)
pg <- ggplot_build(p)

expression_table <- pg$plot$data


ID_list <- c("BOS11A", "BOS12D", "BOS13E" ,"C48A" , "C50C","C51E","C52E","C56E", "ID6","MTS12D","MTS21A","MTS27E","MTS29E","MTS48B","MTS7D", 'CF19', 'CF20', 'CF27', 'CF35', 'CF40')
a=c()
b=c()
for (i in 1:length(ID_list)) {
  ID <- expression_table[expression_table$f %in% ID_list[i],]
  b[length(b) + 1] <- mean(ID$value)
  a[length(a) + 1] <- ID_list[i]
}

mean_table <- data.frame(a,b)

pdf('nCount_SCT_FB_sub_ep_vlnplot.pdf',8,8)
print(VlnPlot(FB_ep_sub, features = 'nCount_SCT', pt.size = 0)+ stat_summary(fun = "mean", geom = "point", color = "red", size = 3))
dev.off()

write.csv(mean_table, 'mean_table_nCount_SCT_sub_ep_vlnplot.csv')

# LTBP2
FB_ep_sub@active.ident <- factor(FB_ep_sub$ID)
p <- VlnPlot2(FB_ep_sub, features = 'LTBP2', pt.size = 0)
pg <- ggplot_build(p)

expression_table <- pg$plot$data


ID_list <- c("BOS11A", "BOS12D", "BOS13E" ,"C48A" , "C50C","C51E","C52E","C56E", "ID6","MTS12D","MTS21A","MTS27E","MTS29E","MTS48B","MTS7D", 'CF19', 'CF20', 'CF27', 'CF35', 'CF40')
a=c()
b=c()
for (i in 1:length(ID_list)) {
  ID <- expression_table[expression_table$f %in% ID_list[i],]
  b[length(b) + 1] <- mean(ID$value)
  a[length(a) + 1] <- ID_list[i]
  
# Convert FB_sub_ep to h5ad object
}

mean_table <- data.frame(a,b)

pdf('LTBP2_FB_sub_ep_vlnplot.pdf',8,8)
print(VlnPlot(FB_ep_sub, features = 'LTBP2', pt.size = 0)+ stat_summary(fun = "mean", geom = "point", color = "red", size = 3))
dev.off()

write.csv(mean_table, 'mean_table_LTBP2_sub_ep_vlnplot.csv')

###############################################################################
saveRDS(FB_ep_sub, 'FB_ep_sub.rds')
