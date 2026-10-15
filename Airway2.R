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
library(nichenetr)
library(tibble)
library(ggplot2)
library(tidyr)
library(ggforce)
plan(multisession(workers = 20))
options(future.globals.maxSize = 8000 * 1024^3)
color.palette = colorRampPalette(c('#0200ad','#fbfcbd','#ff0000'))(255)
setwd("/media/haider/Files/Spatial_Data")
'%ni%' <- function(x,y)!('%in%'(x,y))
library(ggalluvial)

# isolate markers for Mono/Macs cells in the BOS, HC in Epithelial layer
ep_layer <- readRDS("~/Squidpy_Spatial_analysis/ep_obj.rds")
# BOS Ep
BOS_ep <- subset(ep_layer, Group %in% 'BOS')
BOS_marker_genes <- BOS_ep %>% FindAllMarkers(only.pos = TRUE, recorrect_umi=FALSE)
write.csv(BOS_marker_genes, '~/Squidpy_Spatial_analysis/Airway/BOS_Mono_Macs/BOS_Ep_marker_genes.csv')

# HC Ep 
HC_ep <- subset(ep_layer, Group %in% 'HC')
HC_marker_genes <- HC_ep %>% FindAllMarkers(only.pos = TRUE, recorrect_umi=FALSE)
write.csv(HC_marker_genes, '~/Squidpy_Spatial_analysis/Airway/BOS_HC_Mono_Macs/HC_Ep_marker_genes.csv')


#Heatmaps and Volcano Plots (Basal & Mast cells in Epithelial layer)
# COPD vs HC
mast_copd <- subset(ep_layer, subtype %in% c('Mast')) %>% subset(Group %in% c('HC', 'COPD'))
Idents(mast_copd)<-mast_copd$Group
mast_copd.df <- FindMarkers(mast_copd, ident.1 = 'COPD', ident.2 = 'HC',recorrect_umi=FALSE)
mast_copd.df$gene <- rownames(mast_copd.df)

basal_copd <- subset(ep_layer, subtype %in% c('Basal')) %>% subset(Group %in% c('HC', 'COPD'))
Idents(basal_copd)<-basal_copd$Group
basal_copd.df <- FindMarkers(basal_copd, ident.1 = 'COPD', ident.2 = 'HC',recorrect_umi=FALSE)
basal_copd.df$gene <- rownames(basal_copd.df)

# BOS vs HC
#Mast
mast_bos <- subset(ep_layer, subtype %in% c('Mast')) %>% subset(Group %in% c('HC', 'BOS'))
Idents(mast_bos)<-mast_bos$Group
mast_bos.df <- FindMarkers(mast_bos, ident.1 = 'BOS', ident.2 = 'HC',recorrect_umi=FALSE)
mast_bos.df$gene <- rownames(mast_bos.df)

# Basal
basal_bos <- subset(ep_layer, subtype %in% c('Basal')) %>% subset(Group %in% c('HC', 'BOS'))
Idents(basal_bos)<-basal_bos$Group
basal_bos.df <- FindMarkers(basal_bos, ident.1 = 'BOS', ident.2 = 'HC',recorrect_umi=FALSE)
basal_bos.df$gene <- rownames(basal_bos.df)

#Mono_macs
mono_bos <- subset(ep_layer, subtype %in% c("Mono/Macs")) %>% subset(Group %in% c('HC', 'BOS'))
Idents(mono_bos)<-mono_bos$Group
mono_bos.df  <- FindMarkers(mono_bos , ident.1 = 'BOS', ident.2 = 'HC',recorrect_umi=FALSE)
mono_bos.df$gene <- rownames(mono_bos.df)


# CF vs HC
mast_cf <- subset(ep_layer, subtype %in% c('Mast')) %>% subset(Group %in% c('HC', 'CF'))
Idents(mast_cf)<-mast_cf$Group
mast_cf.df <- FindMarkers(mast_cf, ident.1 = 'CF', ident.2 = 'HC',recorrect_umi=FALSE)
mast_cf.df$gene <- rownames(mast_cf.df)

basal_cf <- subset(ep_layer, subtype %in% c('Basal')) %>% subset(Group %in% c('HC', 'CF'))
Idents(basal_cf)<-basal_cf$Group
basal_cf.df <- FindMarkers(basal_cf, ident.1 = 'CF', ident.2 = 'HC',recorrect_umi=FALSE)
basal_cf.df$gene <- rownames(basal_cf.df)

# HC vs Asthma
mast_asthma <- subset(ep_layer, subtype %in% c('Mast')) %>% subset(Group %in% c('HC', 'Asthma'))
Idents(mast_asthma)<-mast_asthma$Group
mast_asthma.df <- FindMarkers(mast_asthma, ident.1 = 'Asthma', ident.2 = 'HC',recorrect_umi=FALSE)
mast_asthma.df$gene <- rownames(mast_asthma.df)

basal_asthma <- subset(ep_layer, subtype %in% c('Basal')) %>% subset(Group %in% c('HC', 'Asthma'))
Idents(basal_asthma)<-basal_asthma$Group
basal_asthma.df <- FindMarkers(basal_asthma, ident.1 = 'Asthma', ident.2 = 'HC',recorrect_umi=FALSE)
basal_asthma.df$gene <- rownames(basal_asthma.df)

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
    scale_color_manual(values =  c("blue", "green",  "firebrick")) +
    theme_grey(base_size = 20)+
    theme(axis.text.x = element_text(face="bold"), text)
  
  print(plot)
  }
setwd("~/Squidpy_Spatial_analysis/Airway/Airway_Mast_Basal/Volcano_plots_mast_basal_HC_vs_allconds")
pdf('basal_mast_COPDvHC_v2.pdf',15,12)
generate_Vol(mast_copd.df)
generate_Vol(basal_copd.df)
dev.off()
s
pdf('basal_mast_CFvHC_v2.pdf',8,8)
generate_Vol(mast_cf.df)
generate_Vol(basal_cf.df)
dev.off()

pdf('basal_mast_BOSvHC_v2.pdf',8,8)
generate_Vol(mast_bos.df)
generate_Vol(basal_bos.df)
dev.off()

pdf('mono_macs_BOSvHC.pdf',15,12)
generate_Vol(mono_bos.df)
dev.off()

pdf('basal_mast_AsthmavHC_v2.pdf',8,8)
generate_Vol(mast_asthma.df)
generate_Vol(basal_asthma.df)
dev.off()

mast_copd_avg <- AverageExpression(mast_copd, return.seurat = TRUE, assays = 'SCT')

basal_copd_avg <- AverageExpression(basal_copd, return.seurat = TRUE, assays = 'SCT')
basal_copd_avg
marker_filter <- function(df) {
  # Read the data from the file
  top_genes <- df %>%
    head(30)
  x <- top_genes %>% dplyr::filter(avg_log2FC>0)
  y <- top_genes %>% dplyr::filter(avg_log2FC<0)
  gene_list <- append(x$gene,y$gene)
  return (gene_list)
}

marker_filter(basal_copd.df) 

DoHeatmap(basal_copd_avg, features = marker_filter(basal_copd.df), draw.lines = F)+
  scale_fill_gradientn(colours = c('blue', 'red'))

?DoHeatmap
DoHeatmap(mast_copd_avg, features = marker_filter(mast_copd.df))+
  scale_fill_gradientn(colours = c('blue', 'red'))
write.csv(basal_copd.df, 'basal_copd.df.csv')
write.csv(mast_copd.df, 'mast_copd.df.csv')
###############################################################################
#Niche Net Receptor Ligand Analysis
basal_mast <- subset(ep_layer, subtype %in% c('Basal', 'Mast'))

#Set up the niche network
organism <- "human"

if(organism == "human"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds"))
} else if(organism == "mouse"){
  lr_network <- readRDS(url("https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"))
  ligand_target_matrix <- readRDS(url("https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"))
  weighted_networks <- readRDS(url("https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"))
  
}

lr_network <- lr_network %>% distinct(from, to)

# Set Receiver
receiver = "Mast"
expressed_genes_receiver <- get_expressed_genes(receiver, basal_mast, pct = 0.05)

all_receptors <- unique(lr_network$to)  
expressed_receptors <- intersect(all_receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>% filter(to %in% expressed_receptors) %>% pull(from) %>% unique()

# Set the Sender Types
sender_celltypes <- c('Basal')

list_expressed_genes_sender <- sender_celltypes %>% unique() %>% lapply(get_expressed_genes, basal_mast, 0.05)
expressed_genes_sender <- list_expressed_genes_sender %>% unlist() %>% unique()

potential_ligands_focused <- intersect(potential_ligands, expressed_genes_sender) 

# Set genes of interest
# condition_oi <-  "COPD"
# condition_reference <- "HC"
# 
# seurat_obj_receiver <- subset(basal_mast, idents = receiver)
# seurat_obj_receiver@active.ident <- as.factor(seurat_obj_receiver$Group)
# DE_table_receiver <-  FindMarkers(object = seurat_obj_receiver,
#                                   ident.1 = condition_oi, ident.2 = condition_reference,
#                                   ,recorrect_umi=FALSE,
#                                   min.pct = 0.05) %>% rownames_to_column("gene")
gene_markers_mb_copd <- FindAllMarkers(subset(basal_mast, Group %in% 'COPD'),recorrect_umi=FALSE)
receiver_genes <- gene_markers_mb_copd[gene_markers_mb_copd$cluster %in% receiver,]
geneset_oi <- receiver_genes  %>% filter(p_val_adj <= 0.05 & avg_log2FC >= 0.25) %>% pull(gene)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

# Set Background genes
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

# Perform analysis and order based on area under precision recall curve
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))

#LogFC of Sender Celltypes
best_upstream_ligands <- ligand_activities %>% top_n(125, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)

# Infer Target Genes and receptors of top ranked ligands
ligand_activity_used <- ligand_activities[ligand_activities$test_ligand %in% rownames(basal_mast),]

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 100) %>%
  bind_rows() %>% drop_na()

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df[ligand_receptor_links_df$from %in% rownames(basal_mast) & ligand_receptor_links_df$to %in% rownames(basal_mast),] ,
  best_upstream_ligands[best_upstream_ligands %in% rownames(basal_mast)],
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Ligands", x_name = "Receptors",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))

# Dotplot of sender-focused approach
gene_markers_all_Ep_Copd <- FindAllMarkers(subset(ep_layer, Group %in% c('COPD')),recorrect_umi=FALSE)
write.csv(gene_markers_all_Ep_Copd, 'gene_markers_all_Ep_Copd.csv')

basal_top_genes <- gene_markers_all_Ep_Copd[gene_markers_all_Ep_Copd$cluster %in% 'Basal',]
basal_top_genes <- basal_top_genes[order(basal_top_genes$pct.1, decreasing = T),] 
basal_top_genes <- basal_top_genes%>% pull(gene)

p_dotplot <- DotPlot(basal_mast,
                     features = rev(basal_top_genes[1:10]), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot

mast_top_genes <- gene_markers_all_Ep_Copd[gene_markers_all_Ep_Copd$cluster %in% 'Mast',]
mast_top_genes <- mast_top_genes[order(mast_top_genes$pct.1, decreasing = T),] 
mast_top_genes <- mast_top_genes%>% pull(gene)

p_dotplot <- DotPlot(basal_mast,
                     features = rev(mast_top_genes[1:10]), cols = "RdYlBu") + 
  coord_flip() +
  scale_y_discrete(position = "right")

p_dotplot

# Get the marker genes for all conditions and just copd
gene_markers_mb <- FindAllMarkers(basal_mast,recorrect_umi=FALSE)
gene_markers_mb_copd <- FindAllMarkers(subset(basal_mast, Group %in% 'COPD'),recorrect_umi=FALSE)
write.csv(gene_markers_mb, 'gene_markers_mast_basal_allconditions.csv')
write.csv(gene_markers_mb, 'gene_markers_mast_basal_copd.csv')


# Alluvial Plot for COPD and HC

ep_layer <- readRDS("~/Squidpy_Spatial_analysis/ep_obj.rds")
copd_ep <- subset(ep_layer, subset=Group %in% c('COPD'))

# Subset out the disease groups of interest and the celltypes
copd_hc_ep <- subset(ep_layer, subset=Group %in% c('COPD', 'HC'))

copd_hc_ep <- subset(copd_hc_ep, subset=subtype %ni% c("DES+_SMC", "gCap", "FCN3+_EC", 
                                                       'SMC', 'FB', 'APOD+_FB', "PLA2G2A+_FB", "DES+_SMC", "aCap", "AT2_EC_DB", "DMBT1+_AT2", 'LQ', 'cDC', 'Peri', 'LEC', 'Cil'))




table <- table(data.frame$Disease_Group, data.frame$Celltypes)
table[,1] <- table[,1]/sum(table[,1])*100
table[,2] <- table[,2]/sum(table[,2])*100


table_high <- table[table$Var1 %ni% c("DES+_SMC", "gCap", "FCN3+_EC", 
                                 'SMC', 'FB', 'APOD+_FB', "PLA2G2A+_FB", "DES+_SMC", "aCap", "AT2_EC_DB", "DMBT1+_AT2", 'LQ','cDC', 'Peri', 'LEC', 'Cil', 'TkC', 'Mast', 'GB','PMN', 
                                 'AM','AT1', "CXCL+_Immune", 'Sec', 'AT2'),]

table_high$Celltypes <- table_high$Var1

table_low <- table[table$Var1 %ni% c("DES+_SMC", "gCap", "FCN3+_EC", 
                                  'SMC', 'FB', 'APOD+_FB', "PLA2G2A+_FB", "DES+_SMC", "aCap", "AT2_EC_DB", "DMBT1+_AT2", 'LQ','cDC', 'Peri', 'LEC', 'Cil', 'TkC', 'Basal', 'Gob', 
                                  'Mono/Macs', 'ThC','BC', 'Perivasc_macs'),]

table_low$Celltypes <- table_low$Var1





#Change the order of celltypes
table_low$Var1 <- factor(table_low$Var1, levels=c('Mast', 'GB','PMN', 'AT1','AM', "CXCL+_Immune", 'Sec', 'AT2'))
table_low$Percent_Abundance <- table_low$Freq
table_high$Var1 <- factor(table_high$Var1, levels=c('Basal', 'Gob', 'Mono/Macs', 'ThC','BC', 'Perivasc_macs'))
table_high$Percent_Abundance <- table_high$Freq


color_pal <- colorRampPalette(c('red', 'darkgreen', 'chocolate', 'magenta', 'cyan','limegreen','orange4', 'maroon', 'gold', 'deeppink1', 'dimgray', 'firebrick', 'purple'))(13)



plot1 <- alluvial_plot(table_high, 'Celltypes with High Abundances in the Airway Epithelial Layer (COPD vs Control)', color_pal[1:6])
plot2 <- alluvial_plot(table_low, 'Celltypes with Low Abundances in the Airway Epithelial Layer (COPD vs Control)', color_pal[7:14])


# Alluvial Plots COPD vs HC

# Alluvial Plot Function

# Produce Smooth curves
ribbon_line <- function(p1, p2, width = 10, len = 100)
{
  if (width > 50) width <- 50
  if (width < 0) width <- 0
  if (p1 < 0) p1 <- 0
  if (p1 > 1) p1 <- 1
  if (p2 < 0) p2 <- 0
  if (p2 > 1) p2 <- 1
  
  yvals <- c(p1, p1, pnorm(seq(-2.5, 2.5, length.out = len)) * (p2 - p1) + p1, p2, p2)
  xvals <- c(0, seq(width, 100 - width, length.out = len + 2), 100)
  list(x = xvals, y = yvals)
}

# Combine the two lines and use them as co-ordinates in a plot
ribbon_df <- function(uppers, lowers, group, width = 10)
{
  data.frame(x    = ribbon_line(uppers[1], uppers[2], width)$x,
             ymax = ribbon_line(uppers[1], uppers[2], width)$y,
             ymin = ribbon_line(lowers[1], lowers[2], width)$y,
             group = group, stringsAsFactors = FALSE)
}

# Create a group of ribbons along with text labels
multi_ribbons <- function(left_bottom, right_bottom, left_top, right_top,
                          groups, width = 10)
{
  if (length(left_bottom) != length(right_bottom) |
      length(left_bottom) != length(left_top)     |
      length(left_top)    != length(right_top)) 
    stop("Left and right columns different length")
  if (length(groups) != length(left_bottom))
    stop("Group length has to be same length as columns")
  
  d <- lapply(seq_along(groups), function(i) {
    ribbon_df(c(left_top[i], right_top[i]),
              c(left_bottom[i], right_bottom[i]),
              groups[i], width)})
  
  left_cols <- lapply(d, function(x) x[1:2,])
  right_cols <- lapply(d, function(x) x[nrow(x) - 1:0,])
  
  res <- list( left = do.call(rbind, left_cols),
               right = do.call(rbind, right_cols),
               bands = do.call(rbind, d))
  
  text_y <- c((res$left$ymax + res$left$ymin)/2,
              (res$right$ymax + res$right$ymin)/2)
  text_x <- c(rep(width / 2, length(res$left$x)), 
              rep(100 - width/2, length(res$left$x)))
  text_labels <- paste0(round(c(res$left$ymax - res$left$ymin,
                                res$right$ymax - res$right$ymin), 3) * 100, "%")
  res$text <- data.frame(x = text_x, y = text_y, labels = text_labels)
  res
}

alluvial <- function(tab, width = 20)
{
  # tab <- table(yvar, xvar)
  # tab <- tab[,colnames(tab) %ni% c("DES+_SMC", "gCap", "FCN3+_EC", 
  #                                  'SMC', 'FB', 'APOD+_FB', "PLA2G2A+_FB", "DES+_SMC", "aCap", "AT2_EC_DB", "DMBT1+_AT2", 'LQ','cDC', 'Peri', 'LEC', 'Cil', 'TkC', 'Sec')]
  
  x_labs <- rownames(tab)
  y_labs <- colnames(tab)
  left <- tab[1,]/sum(tab[1,])
  
  left <- cumsum(sort(left))
  right <- tab[2,]/sum(tab[2,])
  right <- cumsum(sort(right))
  left_lower <- c(0, left[-length(left)])
  names(left_lower) <- names(left)
  right_lower <- c(0, right[-length(right)])
  names(right_lower) <- names(right)
  right <- right[match(names(left), names(right))]
  right_lower <- right_lower[match(names(left), names(right_lower))]
  df_list <- multi_ribbons(left_lower, right_lower, left, right, 
                           names(left), width = 20)
  ggplot(df_list$bands, aes(x = x, ymin = ymin, ymax = ymax, fill = group)) + 
    geom_ribbon(alpha = 0.5) +
    geom_ribbon(alpha = 1, data = df_list$left) +
    geom_ribbon(alpha = 1, data = df_list$right) +
    geom_text(data = df_list$text, inherit.aes = FALSE, colour = "Black",
              aes(x = x, y = y, label = labels), size = 10) +
    geom_text(data = data.frame(x = c(width / 2, 100 - width /2), y = c(1.05, 1.05),
                                labels = factor(x_labs, levels = x_labs)),
              inherit.aes = FALSE,
              mapping = aes(x = x, y = y, label = labels), size = 12) +
    geom_text(data = data.frame(x = rep(-5, length(y_labs)), 
                                y = unique(df_list$text$y[1:(nrow(df_list$text)/2)]), 
                                labs = unique(df_list$bands$group)),
              mapping = aes(x = x, y = y, colour = labs, label = labs),
              inherit.aes = FALSE, size = 10, hjust = 1) +
    scale_fill_manual(values = color_pal) +
    scale_colour_manual(values = color_pal) +
    coord_cartesian(xlim = c(-15, 101)) +
    theme_void() + theme(legend.position = "none")
}
ep_layer <- readRDS("~/Squidpy_Spatial_analysis/ep_obj.rds")
copd_hc_ep <- subset(ep_layer, subset=Group %in% c('HC', 'COPD'))
copd_hc_ep$Group <- factor(copd_hc_ep$Group, levels=c('HC', 'COPD'))
Celltypes <- copd_hc_ep@meta.data$subtype
Disease_Group <- copd_hc_ep@meta.data$Group
data.frame <- data.frame(Disease_Group, Celltypes)

plot3 <- alluvial(data.frame$Disease_Group, data.frame$Celltypes)
plot3
ggsave2('Alluvial_new.pdf', plot3, width = 100, height = 100, limitsize = F)
table(data.frame$Disease_Group, data.frame$Celltypes)

# DotPlot
final_obj <- readRDS("/media/haider/Files/Spatial_Data/annotated_comb_final.rds")
# Remove low quality clusters 
final_obj <- subset(final_obj, subset=subtype %ni% c('LQ',  "AT2_EC_DB" ))
# Get genelists for all genes
all_genes <- final_obj %>% PrepSCTFindMarkers() %>% FindAllMarkers(only.pos = TRUE)

write.csv(all_genes, 'all_genes.csv')

##---Basal---###

# Set markers

markers <- c('KRT7','SCEL','SFTPC', 'SFTPD','DMBT1','MUC5B', 'SERPINA3','LTF','KRT5','KRT15', 'CD24','TMC5',
             'MMRN1', 'PROX1','KDR', 'CA4','ACKR1', 'PLVAP','RAMP2', 'FCN3','PDGFRB', 'KCNK3','ACTA2','MYH11', 'DES', 'LTBP2', 'FGFR4',
             'PLA2G2A','SFRP2', 'APOD','FBN1','CXCL13', 'CXCL10','KIT', 'MS4A2','S100A12','MMP9','MARCO','CD68',
             'CD300E', 'CLEC4E','MS4A4A', 'LYVE1','S100B','FGF2', 'MS4A1', 'FCMR','SEC11C',
             'MZB1', 'CD4','IL7R', 'GZMK', 'CD8A')
# Run Dotplot functions

final_obj <- RenameIdents(final_obj, 'CXCL+_Immune'='CXCL13+_cells')
final_obj <- RenameIdents(final_obj, 'ThC'='CD4 Tcells')
final_obj <- RenameIdents(final_obj, 'TkC'='CD8 Tcells')
                                                                 
final_obj@active.ident <- factor(final_obj@active.ident, levels=c('CD8 Tcells', 'CD4 Tcells', 'GB', 'BC', 'cDC', 'Perivasc_macs', 'Mono/Macs','AM', 'PMN', 'Mast', "CXCL13+_cells",
                                                                  'APOD+_FB', 'PLA2G2A+_FB','FB','DES+_SMC','SMC', 'Peri',
                                                                  'FCN3+_EC', 'gCap','aCap','LEC', 
                                                                  'Cil', 'Basal', 'Sec', 'Gob','DMBT1+_AT2', 'AT2', 'AT1'))                                                              


pdf('all_sample_dotplot.pdf',20, 10)
DotPlot(final_obj, features = markers,  cols = c("blue", "red")) + RotatedAxis()
dev.off()


