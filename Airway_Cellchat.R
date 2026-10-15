library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(future)
library(anndata)
library(plot1cell)
library(SeuratExtend)
library(CellChat)
plan(multisession(workers = 20))
options(future.globals.maxSize = 8000 * 1024^3)
color.palette = colorRampPalette(c('#0200ad','#fbfcbd','#ff0000'))(255)
setwd("~/Squidpy_Spatial_analysis")
################################################################################
data <- read_h5ad("~/Squidpy_Spatial_analysis/Airway_subclass_obj")
`%ni%` <- Negate(`%in%`)
final_obj <- readRDS("/media/haider/Files/Spatial_Data/annotated_comb_final.rds")

#Layers
ep <- data$obs[data$obs$Layers %in% 'Airway_Ep',]

final_obj$X <- sapply(strsplit(as.character(colnames(final_obj)), "_"), `[`, 2)
ep_obj <- subset(final_obj, subset = X %in% ep$cell_ids)
ep_obj@active.assay <- 'SCT'
ep_obj <- subset(ep_obj, subset= ID %ni% c('MTS43D', 'MTS38D'))

saveRDS(ep_obj, 'ep_obj.rds')
levels(factor(ep_obj$subtype))



# Create Individual CellChat objects for Ctrl
ctrl_lis <- subset(ep_obj, Group %in% c("HC"))
data.inpu1 <- GetAssayData(ctrl_lis, assay='SCT', slot='data')
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

ctrl_chat <- identifyOverExpressedGenes(ctrl_chat, do.DE = FALSE)
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
copd_lis <- subset(ep_obj, subset=Group %in% 'COPD')
data.inpu2 <- GetAssayData(copd_lis, assay='SCT', slot='data')
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

copd_chat <- identifyOverExpressedGenes(copd_chat, do.DE = FALSE)
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
copd_chat
################################################################################
# Merge the both Cellchat objects into one single Cellchat object
object.list <- list(ctrl=ctrl_chat, copd=copd_chat)
saveRDS(object.list, 'copd_hc_spatial.rds')

object.list <- readRDS('copd_hc_spatial.rds')

ctrl_chat <-subsetCellChat(ctrl_chat, idents.use = c('Basal', 'Mast'))
copd_chat <-subsetCellChat(copd_chat, idents.use = c('Basal', 'Mast'))
object.list <- list(ctrl=ctrl_chat, copd=copd_chat)
levels(factor(ctrl_chat@idents))
group.celltype <- c("Basal","Mast") 
group.celltype <- factor(group.celltype, levels = c("Basal","Mast"))
object.list <- lapply(object.list, function(x) {mergeInteractions(x, group.celltype)})
cellchat <- mergeCellChat(object.list, add.names=names(object.list))

par(mfrow = c(1, 2), xpd=TRUE)
# compare all the interactions sending from Basal to Mast
#pdf(file=paste0('basal_to_mast', "_chord_diagram.pdf"), width=10, height=10)
for (i in 1:length(object.list)) {
  netVisual_chord_gene(object.list[[i]], sources.use = 2, targets.use = c(1), lab.cex = 0.5, title.name = paste0("Signaling from Basal - ", names(object.list)[i]))}
#dev.off()

saveRDS(cellchat, 'copd_hc_spatial_cellchat_subset.rds')
saveRDS(object.list, file = "cellchat_object.list_subset.rds")