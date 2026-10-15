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
library(zellkonverter)

xenium_layer <- readRDS("/media/haider/Files/Squidpy_Spatial_analysis/final_obj_annotated2.rds")
xenium_layer@images <- list()
xenium_layer$final_annot <- xenium_layer@active.ident

DefaultAssay(xenium_layer) <- "Xenium"

xenium_layer <- JoinLayers(
  xenium_layer,
  assay = "Xenium"
)

sce <- as.SingleCellExperiment(xenium_layer)

writeH5AD(
  sce,
  file = "Whole_Tissue_All_Layers_Lung.h5ad"
)

# Load Seurat object
xenium_layer <- readRDS("/media/haider/Files/Squidpy_Spatial_analysis/final_obj_annotated2.rds")

# Extract counts and set up reticulate objects
anndata <- import("anndata")
scipy <- import("scipy.sparse")
pd <- import("pandas")
np <- import("numpy")

counts <- GetAssayData(xenium_layer, slot = "counts") # or layer = "counts" for Seurat v5
mat <- scipy$csr_matrix(t(as.matrix(counts)))

# Create AnnData
ad <- anndata$AnnData(
  X = mat,
  obs = pd$DataFrame(srt@meta.data),
  var = pd$DataFrame(index = rownames(counts))
)

# Extract spatial coordinates if available
if (length(srt@images) > 0) {
  coords <- srt@images[[1]]@coordinates
  ad$obsm["spatial"] <- np$array(as.matrix(coords[, c("imagerow", "imagecol")]))
}

# Save file
ad$write_h5ad("output_spatial.h5ad")







xenium_layer <- UpdateSeuratObject(xenium_layer)
# DimPlot(xenium_layer, label=T)

xenium_layer@images <- list()
# VolcanoPlots 

# AT1 RAGE Expression COPD VS Non-disease Full tissue

xenium_layer$final_annot <- xenium_layer@active.ident
AT1_copd_hc <- subset(xenium_layer, subset=final_annot %in% 'AT1')
AT1_copd_hc <- subset(AT1_copd_hc, subset=Group %in% c('COPD', 'HC'))

AT1_copd_hc@active.ident <- AT1_copd_hc$Group
AT1_copd_hc.df <- FindMarkers(AT1_copd_hc , ident.1 = 'COPD', ident.2 = 'HC',recorrect_umi=FALSE)
AT1_copd_hc.df$gene <- rownames(AT1_copd_hc.df)
write.csv(AT1_copd_hc.df, "/media/haider/Files/Squidpy_Spatial_analysis/RAGE_ST2_GPR183_Plots/AT1_copd_hc.df.csv")
AT1_copd_hc.df$p_val[AT1_copd_hc.df$p_val==0] <- .Machine$double.xmin



# AT1 RAGE Expression COPD VS Non-disease Full tissue
acap_copd_hc <- subset(xenium_layer, subset=final_annot %in% 'aCap')
acap_copd_hc <- subset(acap_copd_hc, subset=Group %in% c('COPD', 'HC'))

acap_copd_hc@active.ident <- acap_copd_hc$Group
acap_copd_hc.df <- FindMarkers(acap_copd_hc, ident.1 = 'COPD', ident.2 = 'HC',recorrect_umi=FALSE)
acap_copd_hc.df $gene <- rownames(acap_copd_hc.df )
write.csv(acap_copd_hc.df, "/media/haider/Files/Squidpy_Spatial_analysis/RAGE_ST2_GPR183_Plots/acap_copd_hc.df.csv")
acap_copd_hc.df$p_val[acap_copd_hc.df$p_val==0] <- .Machine$double.xmin



# COPD vs CF Mast Gene expression in airway epithelium
mast_copd_cf_airway_ep <- subset(xenium_layer, subset=final_annot %in% 'Mast')
mast_copd_cf_airway_ep   <- subset(mast_copd_cf_airway_ep  , subset=Group %in% c('COPD', 'CF'))
mast_copd_cf_airway_ep   <- subset(mast_copd_cf_airway_ep  , subset= Layers %in% "Airway_Ep")
mast_copd_cf_airway_ep @active.ident <- mast_copd_cf_airway_ep $Group
mast_copd_cf_airway.df <- FindMarkers(mast_copd_cf_airway_ep , ident.1 = 'COPD', ident.2 = 'CF',recorrect_umi=FALSE)
mast_copd_cf_airway.df $gene <- rownames(mast_copd_cf_airway_ep.df  )
# write.csv(mast_copd_cf_airway.df, "/media/haider/Files/Squidpy_Spatial_analysis/RAGE_ST2_GPR183_Plots/mast_copd_cf_airway.df.csv")
# mast_copd_cf_airway.df <- read.csv("/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_V_CF_Volcano_plot/mast_copd_cf_airway_epithelial.csv")
mast_copd_cf_airway_ep.df$p_val[mast_copd_cf_airway_ep.df$p_val==0] <- .Machine$double.xmin


# COPD vs HC Mast Gene expression in airway epithelium
mast_copd_hc_airway_ep <- subset(xenium_layer, subset=final_annot %in% 'Mast')
mast_copd_hc_airway_ep  <- subset(mast_copd_hc_airway_ep , subset=Group %in% c('COPD', 'HC'))
mast_copd_hc_airway_ep  <- subset(mast_copd_hc_airway_ep , subset= Layers %in% "Airway_Ep")
mast_copd_hc_airway_ep@active.ident <- mast_copd_hc_airway_ep $Group
mast_copd_hc_airway_ep.df  <- FindMarkers(mast_copd_hc_airway_ep , ident.1 = 'COPD', ident.2 = 'HC',recorrect_umi=FALSE)
mast_copd_hc_airway_ep.df$gene <- rownames(mast_copd_hc_airway_ep.df )

# # COPD vs Asthma Mast Gene expression in whole
# mast_copd_asthma_airway_whole <- subset(xenium_layer, subset=final_annot %in% 'Mast')
# mast_copd_asthma_airway_whole  <- subset(mast_copd_asthma_airway_whole , subset=Group %in% c('COPD', 'Asthma'))
# mast_copd_asthma_airway_whole@active.ident <- mast_copd_asthma_airway_whole$Group
# mast_copd_asthma_airway_whole.df  <- FindMarkers(mast_copd_asthma_airway_whole , ident.1 = 'Asthma', ident.2 = 'COPD',recorrect_umi=FALSE)
# mast_copd_asthma_airway_whole.df$gene <- rownames(mast_copd_asthma_airway_whole.df )
# mast_copd_asthma_airway_whole.df$p_val[mast_copd_asthma_airway_whole.df$p_val==0] <- .Machine$double.xmin
# 
# 
# # COPD vs Asthma Mast Gene expression in airway epithelium
# mast_copd_asthma_airway_ep <- subset(xenium_layer, subset=final_annot %in% 'Mast')
# mast_copd_asthma_airway_ep  <- subset(mast_copd_asthma_airway_ep , subset=Group %in% c('COPD', 'Asthma'))
# mast_copd_asthma_airway_ep  <- subset(mast_copd_asthma_airway_ep , subset= Layers %in% "Airway_Ep")
# mast_copd_asthma_airway_ep@active.ident <- mast_copd_asthma_airway_ep $Group
# mast_copd_asthma_airway_ep.df  <- FindMarkers(mast_copd_asthma_airway_ep , ident.1 = 'Asthma', ident.2 = 'COPD',recorrect_umi=FALSE)
# mast_copd_asthma_airway_ep.df$gene <- rownames(mast_copd_asthma_airway_ep.df )
# mast_copd_asthma_airway_ep.df$p_val[mast_copd_asthma_airway_ep.df$p_val==0] <- .Machine$double.xmin
# 
# 
# # COPD vs Asthma Mast Gene expression in airway 
# mast_copd_asthma_airway <- subset(xenium_layer, subset=final_annot %in% 'Mast')
# mast_copd_asthma_airway  <- subset(mast_copd_asthma_airway , subset=Group %in% c('COPD', 'Asthma'))
# mast_copd_asthma_airway  <- subset(mast_copd_asthma_airway , subset= Layers %ni% "Non-Airway")
# mast_copd_asthma_airway@active.ident <- mast_copd_asthma_airway $Group
# mast_copd_asthma_airway.df  <- FindMarkers(mast_copd_asthma_airway , ident.1 = 'Asthma', ident.2 = 'COPD',recorrect_umi=FALSE)
# mast_copd_asthma_airway.df$gene <- rownames(mast_copd_asthma_airway.df )
# mast_copd_asthma_airway.df$p_val[mast_copd_asthma_airway.df$p_val==0] <- .Machine$double.xmin


# BC CF vs HC  Gene expression in airway epithelium
BC_cf_hc_airway <- subset(xenium_layer, subset=final_annot %in% 'BC')
BC_cf_hc_airway  <- subset(BC_cf_hc_airway , subset=Group %in% c('CF', 'HC'))
BC_cf_hc_airway  <- subset(BC_cf_hc_airway , subset= Layers %in% "Non-Airway")
BC_cf_hc_airway@active.ident <- BC_cf_hc_airway $Group
BC_cf_hc_airway.df  <- FindMarkers(BC_cf_hc_airway , ident.1 = 'CF', ident.2 = 'HC',recorrect_umi=FALSE)
BC_cf_hc_airway.df$gene <- rownames(BC_cf_hc_airway.df )

# Mono BOS vs HC  Gene expression in airway epithelium
mono_bos_hc_airway <- subset(xenium_layer, subset=final_annot %in% c('Mono/Macs', 'AM',"Perivasc_macs"))
mono_bos_hc_airway   <- subset(mono_bos_hc_airway  , subset=Group %in% c('BOS', 'HC'))
# mono_bos_hc_airway   <- subset(mono_bos_hc_airway  , subset= Layers %ni%  "Non-Airway")
mono_bos_hc_airway @active.ident <- mono_bos_hc_airway  $Group
mono_bos_hc_airway.df  <- FindMarkers(mono_bos_hc_airway  , ident.1 = 'BOS', ident.2 = 'HC',recorrect_umi=FALSE)
mono_bos_hc_airway.df$gene <- rownames(mono_bos_hc_airway.df )
mono_bos_hc_airway.df$p_val[mono_bos_hc_airway.df$p_val==0] <- .Machine$double.xmin

# Mast Pseudo bulk analysis COPD vs Asthma
mast_copd_asthma_airway_pseudo <- subset(xenium_layer, subset = final_annot %in% 'Mast')
mast_copd_asthma_airway_pseudo <- subset(mast_copd_asthma_airway_pseudo, subset = Group %in% c('COPD', 'Asthma'))
mast_copd_asthma_airway_pseudo <- subset(mast_copd_asthma_airway_pseudo, subset = Layers %ni% c('Non-Airway'))

mast_copd_asthma_airway_pseudo <- AggregateExpression(
  mast_copd_asthma_airway_pseudo,
  assays = "Xenium",
  return.seurat = TRUE,
  group.by = c("Group", "ID", "final_annot")   # <-- ID is your real patient/specimen column
)

mast_copd_asthma_airway_pseudo$celltype.stim <- paste(mast_copd_asthma_airway_pseudo$final_annot, mast_copd_asthma_airway_pseudo$Group, sep = "_")
Idents(mast_copd_asthma_airway_pseudo) <- "celltype.stim"
# Deseq2 bulk analysis

mast_copd_asthma_airway_pseudo.df <- FindMarkers(object = mast_copd_asthma_airway_pseudo, 
                            ident.1 = "Mast_Asthma", 
                            ident.2 = "Mast_COPD",
                            test.use = "DESeq2")
head(mast_copd_asthma_airway_pseudo.df, n = 15)


# Mast cells airway vs non airway healthy
xenium_layer$Layers_2 <- xenium_layer$Layers
xenium_layer$Layers_2[xenium_layer$Layers_2 %in% c('Airway_Ep', 'Sub_Ep', 'SMG')] <- paste('Airway')

xenium_layer$Group <- as.character(xenium_layer$Group)
xenium_layer$Group[xenium_layer$Group %in% c('HC')] <- paste0('ND')
xenium_layer$Group <- as.factor(xenium_layer$Group)

mast_airway_par_ctrl <- subset(xenium_layer, subset=final_annot %in% ('Mast')) %>% subset(subset=Group %in% ('ND'))
mast_airway_par_ctrl@active.ident <- factor(mast_airway_par_ctrl$Layers_2)
mast_airway_par_ctrl.df  <- FindMarkers(mast_airway_par_ctrl , ident.1 = 'Airway', ident.2 = 'Non-Airway',recorrect_umi=FALSE)
mast_airway_par_ctrl.df$gene <- rownames(mast_airway_par_ctrl.df )
mast_airway_par_ctrl.df$p_val_adj <- p.adjust(mast_airway_par_ctrl.df$p_val, method = "BH")

# Mast cells airway vs non airway copd
mast_airway_par_copd <- subset(xenium_layer, subset=final_annot %in% ('Mast')) %>% subset(subset=Group %in% ('COPD'))
mast_airway_par_copd@active.ident <- factor(mast_airway_par_copd $Layers_2)
mast_airway_par_copd.df  <- FindMarkers(mast_airway_par_copd , ident.1 = 'Airway', ident.2 = 'Non-Airway',recorrect_umi=FALSE)
mast_airway_par_copd.df$gene <- rownames(mast_airway_par_copd.df)
mast_airway_par_copd.df$p_val_adj <- p.adjust(mast_airway_par_copd.df$p_val, method = "BH")


generate_Vol <- function(res.df){
  
  res.df$Regulation <- "NS"
  
  res.df$Regulation[res.df$avg_log2FC > 1 & res.df$p_val_adj < 0.05] <- "UP"
  res.df$Regulation[res.df$avg_log2FC < -1 & res.df$p_val_adj < 0.05] <- "Down"
  
  up_genes   <- dplyr::filter(res.df, Regulation == "UP")
  down_genes <- dplyr::filter(res.df, Regulation == "Down")
  
  print(up_genes)
  
  plot <- ggplot(res.df, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
    
    geom_point(aes(color = Regulation), alpha = 0.4, size = 2) +
    
    geom_vline(xintercept = c(-1 ,1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    
    geom_point(
      data = up_genes,
      aes(color = Regulation),
      shape = 21,
      size = 3,
      fill = "firebrick"
    ) +
    
    geom_point(
      data = down_genes,
      aes(color = Regulation),
      shape = 21,
      size = 3,
      fill = "blue"
    ) +
    # ggrepel::geom_text_repel(
    #   data = up_genes[up_genes$gene %in% 'IL13',],
    #   aes(label = gene),
    #   size = 2
    # ) +
    # ggrepel::geom_text_repel(
    #   data = up_genes,
    #   aes(label = gene),
    #   size = 2
    # ) +
    # 
    # 
    # ggrepel::geom_text_repel(
    #   data = down_genes,
    #   aes(label = gene),
    #   size = 2
    # ) +


    scale_color_manual(values = c(
      "UP" = "firebrick",
      "Down" = "blue",
      "NS" = "grey"
    )) +
    
    theme_bw(base_size = 8)+
    theme(
      panel.grid = element_blank(),
      legend.position = c(0.92,0.85)
    )
  print(plot)
}



pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_V_CF_Volcano_plot/mast_copd_cf_airway_epithelial_volcano_plot.pdf", 10, 6)
generate_Vol(mast_copd_cf_airway.df)
dev.off()

write.csv(mast_copd_cf_airway.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_V_CF_Volcano_plot/mast_copd_cf_airway_epithelial.csv")



pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_HC/mast_copd_cf_airway_epithelial_volcano_plot.pdf", 10, 6)
generate_Vol(mast_copd_hc_airway.df)
dev.off()

pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_HC/mast_copd_cf_airway_epithelial_volcano_plot.pdf", 10, 6)
generate_Vol(mast_copd_hc_airway.df)
dev.off()


write.csv(mast_copd_hc_airway.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_HC/mast_airway_hc_copd.df.csv")



generate_Vol(BC_cf_hc_airway.df)


# whole tissue Mono BOS vs HC
all_macs_bos_hc_whole.df <- read.csv("/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/AM_Perivasc_MonoMac_bos_hc_whole_tissue.csv")
all_macs_bos_hc_whole.df$p_val_adj <- p.adjust(all_macs_bos_hc_whole.df$p_val, method = "BH")
write.csv(all_macs_bos_hc_whole.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/AM_Perivasc_MonoMac_bos_hc_whole_tissue+bh_fdr.csv")

all_macs_bos_hc_whole.df$p_val_adj[all_macs_bos_hc_whole.df$p_val_adj==0] <- .Machine$double.xmin


pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/AM_Perivasc_MonoMac_bos_hc_whole_tissue_nolabel_fdr.pdf", 6,4)
generate_Vol(all_macs_bos_hc_whole.df) 
dev.off()

pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/AM_Perivasc_MonoMac_bos_hc_whole_tissue_fdr.pdf", 6,4)
generate_Vol(all_macs_bos_hc_whole.df) 
dev.off()

write.csv(all_macs_bos_hc_whole.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/AM_Perivasc_MonoMac_bos_hc_whole_tissue.csv")


# airway Mono BOS vs HC
mono_bos_hc_airway.df <- read.csv("/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/mono_mac_bos_hc_airway.csv")
mono_bos_hc_airway.df$p_val_adj <- p.adjust(mono_bos_hc_airway.df$p_val, method = "BH")
write.csv(mono_bos_hc_airway.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/mono_mac_bos_hc_airway_bh_fdr.csv")


mono_bos_hc_airway.df $p_val_adj[mono_bos_hc_airway.df $p_val_adj==0] <- .Machine$double.xmin


pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/mono_mac_bos_hc_airway_nolabel_fdr.pdf", 6,4)
generate_Vol(mono_bos_hc_airway.df) 
dev.off()

pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/mono_mac_bos_hc_airway_fdr.pdf", 6,4)
generate_Vol(mono_bos_hc_airway.df) 
dev.off()

write.csv(mono_bos_hc_airway.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/mono_mac_bos_hc_airway.csv")

# whole Mono BOS vs HC

mono_bos_hc_whole.df <- read.csv("/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/mono_mac_bos_hc_whole_tissue.csv")
mono_bos_hc_whole.df$p_val_adj <- p.adjust(mono_bos_hc_whole.df$p_val, method = "BH")
write.csv(mono_bos_hc_whole.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/mono_mac_bos_hc_whole_tissue_bh_fdr.csv")


mono_bos_hc_whole.df  $p_val_adj[mono_bos_hc_whole.df  $p_val_adj==0] <- .Machine$double.xmin


pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/mono_mac_bos_hc_whole_nolabel_fdr.pdf", 6,4)
generate_Vol(mono_bos_hc_whole.df ) 
dev.off()

pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/mono_mac_bos_hc_whole_fdr.pdf", 6,4)
generate_Vol(mono_bos_hc_whole.df ) 
dev.off()

write.csv(mono_bos_hc_whole.df , "/media/haider/Files/Squidpy_Spatial_analysis/Mono_Mac_BOS_vs_HC/mono_mac_bos_hc_airway.csv")

##################################################################################

BC_cf_hc_whole.df <- read.csv("/media/haider/Files/Squidpy_Spatial_analysis/BC_CF_vs_HC/bc_cf_hc_whole_tissue.csv")
BC_cf_hc_whole.df$p_val_adj <- p.adjust(BC_cf_hc_whole.df$p_val, method = "BH")
write.csv(BC_cf_hc_whole.df, "/media/haider/Files/Squidpy_Spatial_analysis/BC_CF_vs_HC/bc_cf_hc_whole_tissue_bh_fdr.csv")

# whole tissue BC CF vs HC
pdf("/media/haider/Files/Squidpy_Spatial_analysis/BC_CF_vs_HC/bc_cf_hc_whole_tissue.pdf_cutoff_1_logfc_fdr", 6,4)
generate_Vol(BC_cf_hc_whole.df ) + xlim(-4,4)
dev.off()

pdf("/media/haider/Files/Squidpy_Spatial_analysis/BC_CF_vs_HC/bc_cf_hc_whole_tissue_nolabel_cutoff_1logfc_fdr.pdf", 6,4)
generate_Vol(BC_cf_hc_whole.df ) + xlim(-4,4)
dev.off()

write.csv(BC_cf_hc_whole.df, "/media/haider/Files/Squidpy_Spatial_analysis/BC_CF_vs_HC/bc_cf_hc_whole_tissue.csv")


# airway  BC CF vs HC
BC_cf_hc_airway.df <- read.csv("/media/haider/Files/Squidpy_Spatial_analysis/BC_CF_vs_HC/bc_cf_hc_airway.csv")
BC_cf_hc_airway.df$p_val_adj <- p.adjust(BC_cf_hc_airway.df$p_val, method = "BH")
write.csv(BC_cf_hc_airway.df, "/media/haider/Files/Squidpy_Spatial_analysis/BC_CF_vs_HC/bc_cf_hc_airway_bh_fdr.csv")

pdf("/media/haider/Files/Squidpy_Spatial_analysis/BC_CF_vs_HC/bc_cf_hc_airway_logfc_1_fdr.pdf", 6,4)
generate_Vol(BC_cf_hc_airway.df) 
dev.off()

pdf("/media/haider/Files/Squidpy_Spatial_analysis/BC_CF_vs_HC/bc_cf_hc_airway_nolabel_logfc_1_fdr.pdf", 6,4)
generate_Vol(BC_cf_hc_airway.df)  
dev.off()

write.csv(BC_cf_hc_airway.df, "/media/haider/Files/Squidpy_Spatial_analysis/BC_CF_vs_HC/bc_cf_hc_airway.csv")

#Airway Ep Mast HC vs COPD
mast_copd_hc_airway_ep.df <- read.csv("/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_HC/mast_airway_ep_hc_copd.df.csv")
mast_copd_hc_airway_ep.df$p_val_adj <- p.adjust(mast_copd_hc_airway_ep.df$p_val, method = "BH")
write.csv(mast_copd_hc_airway_ep.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_HC/mast_airway_ep_hc_copd.df_bh_fdr.csv")

pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_HC/mast_copd_hc_airway_epithelial_volcano_plot_logfc1_label_fdr.pdf", 6, 4)
generate_Vol(mast_copd_hc_airway_ep.df)
dev.off()

pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_HC/mast_copd_hc_airway_epithelial_volcano_plot_logfc1_nolabel_fdr.pdf", 6, 4)
generate_Vol(mast_copd_hc_airway_ep.df)
dev.off()

write.csv(mast_copd_hc_airway_ep.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_HC/mast_airway_ep_hc_copd.df.csv")



# Airway Ep Mast CF vs COPD
mast_copd_cf_airway.df <- read.csv("/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_V_CF_Volcano_plot/mast_copd_cf_airway_epithelial.csv")
mast_copd_cf_airway.df$p_val_adj <- p.adjust(mast_copd_cf_airway.df$p_val, method = "BH")
write.csv(mast_copd_cf_airway.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_V_CF_Volcano_plot/mast_copd_cf_airway_epithelial_bh_fdr.csv")

mast_copd_cf_airway.df$p_val_adj[mast_copd_cf_airway.df$p_val_adj==0] <- .Machine$double.xmin


pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_V_CF_Volcano_plot/mast_copd_cf_airway_epithelial_label_fdr.pdf", 6, 4)
generate_Vol(mast_copd_cf_airway.df)
dev.off()

pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_COPD_V_CF_Volcano_plot/mast_copd_cf_airway_epithelial_nolabel_fdr.pdf", 6, 4)
generate_Vol(mast_copd_cf_airway.df)
dev.off()


######################################################################################################################################
# Airway Mast COPD vs Asthma
pdf("/media/haider/Files/Squidpy_Spatial_analysis/Volcano_Mast_COPD_Asthma/copd_asthma_mast_airway_nolabel.pdf", 12, 6)
generate_Vol(mast_copd_asthma_airway.df)
dev.off()

# Airway Ep Mast COPD vs Asthma
pdf("/media/haider/Files/Squidpy_Spatial_analysis/Volcano_Mast_COPD_Asthma/copd_asthma_mast_airway_ep_label.pdf",  12, 6)
generate_Vol(mast_copd_asthma_airway_ep.df)
dev.off()

# Whole Mast COPD vs Asthma
pdf("/media/haider/Files/Squidpy_Spatial_analysis/Volcano_Mast_COPD_Asthma/copd_asthma_mast_airway_whole_nolabel.pdf",10, 6)
generate_Vol(mast_copd_asthma_airway_whole.df)
dev.off()

write.csv(mast_copd_asthma_airway_ep.df, "/media/haider/Files/Squidpy_Spatial_analysis/Volcano_Mast_COPD_Asthma/copd_asthma_mast_airway_ep.csv")
write.csv(mast_copd_asthma_airway.df, "/media/haider/Files/Squidpy_Spatial_analysis/Volcano_Mast_COPD_Asthma/copd_asthma_mast_airway.csv")
write.csv(mast_copd_asthma_airway_whole.df, "/media/haider/Files/Squidpy_Spatial_analysis/Volcano_Mast_COPD_Asthma/copd_asthma_mast_whole.csv")

mast_copd_asthma_airway_pseudo.df$gene <- rownames(mast_copd_asthma_airway_pseudo.df)
generate_Vol(mast_copd_asthma_airway_pseudo.df)

################################################################################
# GPR183 Mast cells 
#COPD
pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_Non_airway_vs_Airway_COPD_ND/mast_airway_par_copd_volcano_plot_label.pdf", 6, 4)
generate_Vol(mast_airway_par_copd.df)
dev.off()

pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_Non_airway_vs_Airway_COPD_ND/mast_airway_par_copd_volcano_plot_nolabel.pdf", 6, 4)
generate_Vol(mast_airway_par_copd.df)
dev.off()
write.csv(mast_airway_par_copd.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mast_Non_airway_vs_Airway_COPD_ND/mast_airway_par_copd.csv")

#Healthy
pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_Non_airway_vs_Airway_COPD_ND/mast_airway_par_healthy_volcano_plot_label.pdf", 6, 4)
generate_Vol(mast_airway_par_ctrl.df)
dev.off()

pdf("/media/haider/Files/Squidpy_Spatial_analysis/Mast_Non_airway_vs_Airway_COPD_ND/mast_airway_par_healthy_volcano_plot_nolabel.pdf", 6, 4)
generate_Vol(mast_airway_par_ctrl.df)
dev.off()

write.csv(mast_airway_par_ctrl.df, "/media/haider/Files/Squidpy_Spatial_analysis/Mast_Non_airway_vs_Airway_COPD_ND/mast_airway_par_healthy.csv")
