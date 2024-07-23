#############################################################################
## Weller, Bartok, McGinnis et al., 2024 companion analysis code ############
## Translation dysregulation in cancer as a source for targetable antigens ##
## July 23rd, 2024; Chris McGinnis, PhD; Stanford Univeristy, Satpathy Lab ##
#############################################################################

library(Seurat)
library(ggplot2)
library(ComplexHeatmap)
library(reshape2)
library(viridis)
library(speckle)
library(RColorBrewer)
library(CellChat)
library(ComplexHeatmap)

#######################
## Main Text Figures ##
#######################

load("tumor_vol_summary.Robj")
load("seu_ys_final.Robj")
load('anno_freq_all.Robj')
load('seu_ys_cd8t.Robj')
load('umap_cd8t.Robj')
load('anno_freq_cd8t.Robj')
load('exh_zscores.Robj')
load('codex_summary.Robj')
load('codex_data.Robj')
load('cc_obj_list.Robj')

## Fig. 6A: Tumor growth curves
ggplot(tumor_vol_summary, aes(x=time, y=mean, color=clone)) +
  geom_point() + geom_line() + theme_classic() + scale_color_manual(values=c('#ff8989','#DC0000FF','#8C4484')) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.5, color='black') +
  theme(legend.position = 'none')

## Fig. 6C: Immune cell type UMAP
DimPlot(seu_ys_final, group.by = 'celltype_fig2', cols=c('lightcoral','red','darkred','navy','lightsalmon','steelblue3','maroon','grey60','pink3'))+NoLegend()+NoAxes()+theme(plot.title = element_blank())

## Fig. 6D: CD8T proportion timepoint barchart
ggplot(anno_freq_all[grep('cd8', anno_freq_all$celltype), ], aes(x=day, fill=sample, y=freq)) + geom_col(width = 0.8, color='black', position=position_dodge()) + theme_classic() + scale_fill_manual(values=alpha(c('#8C4484','#ff8989','#DC0000FF','#8C4484','#ff8989','#DC0000FF'),0.8)) +
  theme(legend.position = 'none', axis.title.x = element_blank(), axis.title.y = element_blank()) 

## Fig. 6E: CD8T subtype UMAP
DimPlot(seu_ys_cd8t, group.by = 'subtype', cols=c('red','lightcoral','black','grey40','darkred','lightsalmon')) + NoLegend() + NoAxes() + theme(plot.title = element_blank())

## Fig. 6F: CD8T subtype annotation dotplot
anno_markers_cd8 <- c('Ccr7','Sell','Satb1','Il7r','Tcf7','Cxcr3','Gzmk','Gzmb','Ccl5','Tox','Lag3','Pdcd1','Xcl1','Tnfrsf9','Rel','Mki67','Hells','Pclaf') 
g <- DotPlot(seu_ys_cd8t, group.by = 'subtype', features=anno_markers_cd8, cols='RdBu', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
g$data$id <- factor(g$data$id, levels=rev(c('cd8_naive','cd8_memlike','cd8_efflike','cd8_tex','cd8_tpex','cd8_prolif')))
print(g)

## Fig. 6G: CD8T sample UMAPs
ggplot(umap_cd8t, aes(x=umap_1,y=umap_2)) + geom_point(data=umap_cd8t, size=3.95, color="black") + geom_point(data=umap_cd8t, size=3.55, color="gray90") +
  geom_point(data = umap_cd8t[which(umap_cd8t$sample =='wt19_early'),], color=alpha('#8C4484',0.8)) + theme_void() 
ggplot(umap_cd8t, aes(x=umap_1,y=umap_2)) + geom_point(data=umap_cd8t, size=3.95, color="black") + geom_point(data=umap_cd8t, size=3.55, color="gray90") +
  geom_point(data = umap_cd8t[which(umap_cd8t$sample =='ko13_early'),], color=alpha('#ff8989',0.8)) + theme_void() 
ggplot(umap_cd8t, aes(x=umap_1,y=umap_2)) + geom_point(data=umap_cd8t, size=3.95, color="black") + geom_point(data=umap_cd8t, size=3.55, color="gray90") +
  geom_point(data = umap_cd8t[which(umap_cd8t$sample =='ko18_early'),], color=alpha('#DC0000FF',0.8)) + theme_void() 
ggplot(umap_cd8t, aes(x=umap_1,y=umap_2)) + geom_point(data=umap_cd8t, size=3.95, color="black") + geom_point(data=umap_cd8t, size=3.55, color="gray90") +
  geom_point(data = umap_cd8t[which(umap_cd8t$sample =='wt19_late'),], color=alpha('#8C4484',0.8)) + theme_void() 
ggplot(umap_cd8t, aes(x=umap_1,y=umap_2)) + geom_point(data=umap_cd8t, size=3.95, color="black") + geom_point(data=umap_cd8t, size=3.55, color="gray90") +
  geom_point(data = umap_cd8t[which(umap_cd8t$sample =='ko13_late'),], color=alpha('#ff8989',0.8)) + theme_void() 
ggplot(umap_cd8t, aes(x=umap_1,y=umap_2)) + geom_point(data=umap_cd8t, size=3.95, color="black") + geom_point(data=umap_cd8t, size=3.55, color="gray90") +
  geom_point(data = umap_cd8t[which(umap_cd8t$sample =='ko18_late'),], color=alpha('#DC0000FF',0.8)) + theme_void() 

## Fig. 6H: CD8 subtype proportion timepoint barchart
ggplot(anno_freq_cd8t[grep('tex', anno_freq_cd8t$celltype), ], aes(x=day, fill=sample, y=freq)) + geom_col(width = 0.8, color='black', position=position_dodge()) + theme_classic() + scale_fill_manual(values=alpha(c('#8C4484','#ff8989','#DC0000FF','#8C4484','#ff8989','#DC0000FF'),0.8)) +
  theme(legend.position = 'none', axis.title.x = element_blank(), axis.title.y = element_blank()) 
ggplot(anno_freq_cd8t[grep('act', anno_freq_cd8t$celltype), ], aes(x=day, fill=sample, y=freq)) + geom_col(width = 0.8, color='black', position=position_dodge()) + theme_classic() + scale_fill_manual(values=alpha(c('#8C4484','#ff8989','#DC0000FF','#8C4484','#ff8989','#DC0000FF'),0.8)) +
  theme(legend.position = 'none', axis.title.x = element_blank(), axis.title.y = element_blank()) 

## Fig. 6I: CD8 exhaustion marker z-score heatmap
Heatmap(exh_zscores, cluster_columns = T, cluster_rows = F, col=rev(brewer.pal(name='RdBu',n=9)))

## Fig. 6J: CODEX analysis
ggplot(codex_summary[which(codex_summary$measure == 'cd8'),], aes(x=geno, fill=geno, y=mean)) + geom_col(color='black', width=0.5) + theme_classic() + scale_fill_manual(values=alpha(c('#8C4484','#DC0000FF'),0.8)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) + theme(legend.position = 'none', axis.line = element_line(size=1), axis.ticks = element_line(size=1), axis.ticks.length=unit(.2, "cm")) + 
  geom_point(data=codex_data[which(codex_data$measure == 'cd8'),], aes(x=geno, y=value)) 

ggplot(codex_summary[which(codex_summary$measure == 'lag3'),], aes(x=geno, fill=geno, y=mean)) + geom_col(color='black', width=0.5) + theme_classic() + scale_fill_manual(values=alpha(c('#8C4484','#DC0000FF'),0.8)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) + theme(legend.position = 'none', axis.line = element_line(size=1), axis.ticks = element_line(size=1), axis.ticks.length=unit(.2, "cm")) + 
  geom_point(data=codex_data[which(codex_data$measure == 'lag3'),], aes(x=geno, y=value)) 

ggplot(codex_summary[which(codex_summary$measure == 'ki67'),], aes(x=geno, fill=geno, y=mean)) + geom_col(color='black', width=0.5) + theme_classic() + scale_fill_manual(values=alpha(c('#8C4484','#DC0000FF'),0.8)) +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=0.2) + theme(legend.position = 'none', axis.line = element_line(size=1), axis.ticks = element_line(size=1), axis.ticks.length=unit(.2, "cm")) + 
  geom_point(data=codex_data[which(codex_data$measure == 'ki67'),], aes(x=geno, y=value)) 

t.test(x=codex_data$value[which(codex_data$geno == 'WT' & codex_data$measure == 'cd8')],
       y=codex_data$value[which(codex_data$geno == 'KO' & codex_data$measure == 'cd8')])
t.test(x=codex_data$value[which(codex_data$geno == 'WT' & codex_data$measure == 'lag3')],
       y=codex_data$value[which(codex_data$geno == 'KO' & codex_data$measure == 'lag3')])
t.test(x=codex_data$value[which(codex_data$geno == 'WT' & codex_data$measure == 'ki67')],
       y=codex_data$value[which(codex_data$geno == 'KO' & codex_data$measure == 'ki67')])

## Fig. 6K: CellChat IFNG signaling network
pathways.show <- 'IFN-II'
weight.max <- getMaxWeight(cc_obj_list, slot.name = c("netP"), attribute = pathways.show) 
netVisual_individual(cc_obj_list[[1]], edge.weight.max = weight.max[1], signaling = pathways.show, pairLR.use = 'IFNG_IFNGR1_IFNGR2', color.use=c('navy','darkred','steelblue3','maroon','red','lightsalmon','lightcoral','pink3'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01, )
netVisual_individual(cc_obj_list[[3]], edge.weight.max = weight.max[1], signaling = pathways.show, pairLR.use = 'IFNG_IFNGR1_IFNGR2', color.use=c('navy','darkred','steelblue3','maroon','red','lightsalmon','lightcoral','pink3'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01)
netVisual_individual(cc_obj_list[[5]], edge.weight.max = weight.max[1], signaling = pathways.show, pairLR.use = 'IFNG_IFNGR1_IFNGR2', color.use=c('navy','darkred','steelblue3','maroon','red','lightsalmon','lightcoral','pink3'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01)
netVisual_individual(cc_obj_list[[2]], edge.weight.max = weight.max[1], signaling = pathways.show, pairLR.use = 'IFNG_IFNGR1_IFNGR2', color.use=c('navy','darkred','steelblue3','maroon','red','lightsalmon','lightcoral','pink3'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01)
netVisual_individual(cc_obj_list[[4]], edge.weight.max = weight.max[1], signaling = pathways.show, pairLR.use = 'IFNG_IFNGR1_IFNGR2', color.use=c('navy','darkred','steelblue3','maroon','red','lightsalmon','lightcoral','pink3'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01)
netVisual_individual(cc_obj_list[[6]], edge.weight.max = weight.max[1], signaling = pathways.show, pairLR.use = 'IFNG_IFNGR1_IFNGR2', color.use=c('navy','darkred','steelblue3','maroon','red','lightsalmon','lightcoral','pink3'), vertex.label.color = 'white', arrow.size=0.2, thresh = 0.01)

##########################
## Supplemental Figures ##
##########################

load('anno_markers_all.Robj')
load('anno_freq_all.Robj')
load('prop_test_all.Robj')
load('anno_freq_cd8t.Robj')
load('prop_test_cd8t.Robj')
load('seu_ys_nk.Robj')
load('anno_freq_nk.Robj')

## Fig. S10A: Immune cell type annotation dotplot
g <- DotPlot(seu_ys_final, group.by = 'celltype_fig2', features=anno_markers_all, cols='RdBu', dot.scale = 3) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank())
g$data$id <- factor(g$data$id, levels=rev(c('b','cd4t','cd8t','ilc','treg','nk','myl','dc','prolif')))
print(g)

## Fig. S10B: Immune cell type proportion timepoint barchart
ggplot(anno_freq_all[grep('early',anno_freq_all$sample), ], aes(x=celltype, fill=sample, y=freq)) + geom_col(color='black', position = position_dodge()) + theme_classic() + scale_fill_manual(values=alpha(c('#8C4484','#ff8989','#DC0000FF'),0.8)) + theme(legend.position = 'none', axis.title.x = element_blank(), axis.title.y = element_blank()) + scale_y_continuous(position = "right", limits = c(0,0.45)) 
ggplot(anno_freq_all[grep('late',anno_freq_all$sample), ], aes(x=celltype, fill=sample, y=freq)) + geom_col(color='black', position = position_dodge()) + theme_classic() + scale_fill_manual(values=alpha(c('#8C4484','#ff8989','#DC0000FF'),0.8)) + theme(legend.position = 'none', axis.title.x = element_blank(), axis.title.y = element_blank()) + scale_y_continuous(position = "right", limits = c(0,0.45)) 

## Fig. S10C: Immune cell type proportion change p-value heatmap
prop_test_all_hm <- acast(prop_test_all, comparison~celltype, value.var="pval_bin")
prop_test_all_hm <- prop_test_all_hm[rev(c('wt_ko13_early','wt_ko18_early','wt_ko13_late','wt_ko18_late')), c('b','cd4t','cd8t','ilc','treg','nk','myl','dc','prolif')]
Heatmap(prop_test_all_hm, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, col = c(alpha('black',0.8),alpha(viridis(4)[4],0.8)))

## Fig. S10D: CD8T subtype proportion timepoint barchart
ggplot(anno_freq_cd8t[grep('early',anno_freq_cd8t$sample), ], aes(x=celltype, fill=sample, y=freq)) + geom_col(color='black', position = position_dodge()) + theme_classic() + scale_fill_manual(values=alpha(c('#8C4484','#ff8989','#DC0000FF'),0.8)) + theme(legend.position = 'none', axis.title.x = element_blank(), axis.title.y = element_blank()) + scale_y_continuous(position = "right", limits = c(0,0.4)) 
ggplot(anno_freq_cd8t[grep('late',anno_freq_cd8t$sample), ], aes(x=celltype, fill=sample, y=freq)) + geom_col(color='black', position = position_dodge()) + theme_classic() + scale_fill_manual(values=alpha(c('#8C4484','#ff8989','#DC0000FF'),0.8)) + theme(legend.position = 'none', axis.title.x = element_blank(), axis.title.y = element_blank()) + scale_y_continuous(position = "right", limits = c(0,0.4)) 

## Fig. S10E: CD8T p-value heatmap
prop_test_cd8t_hm <- acast(prop_test_cd8t, comparison~celltype, value.var="pval_bin")
prop_test_cd8t_hm <- prop_test_cd8t_hm[rev(c('wt_ko13_early','wt_ko18_early','wt_ko13_late','wt_ko18_late')), c('cd8_naive','cd8_act','cd8_gzmk','cd8_tex','cd8_tpex','cd8_prolif')]
Heatmap(prop_test_cd8t_hm, cluster_rows = F, cluster_columns = F, show_heatmap_legend = F, col = rev(c(alpha('black',0.8),alpha(viridis(4)[4],0.8))))

## Fig. S10F: NK subtype UMAP and annotation marker feature plots
DimPlot(seu_ys_nk, group.by = 'subtype2', cols=c('seagreen','palegreen3','grey40'), pt.size = 1) + NoLegend() + NoAxes() + theme(plot.title = element_blank())
FeaturePlot(seu_ys_nk, 'Prf1', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis()+theme(plot.title = element_blank())
FeaturePlot(seu_ys_nk, 'Ctla2a', max.cutoff = 'q95', pt.size = 1) + NoLegend() + NoAxes() + scale_color_viridis()+theme(plot.title = element_blank())

## Fig. S10G: Cytotoxic NK cell subtype proporton and cytotoxicity score timepoint barchart and violin plot
ggplot(anno_freq_nk[grep('nk_cyto', anno_freq_nk$celltype), ], aes(x=clone, fill=day_clone, y=freq)) + geom_col(width = 0.8, color='black', position=position_dodge()) + theme_classic() + scale_fill_manual(values=alpha(c('#ff8989','#DC0000FF','#8C4484','#ff8989','#DC0000FF','#8C4484'),0.8)) +
  theme(legend.position = 'none', axis.title.x = element_blank(), axis.title.y = element_blank())

# cytotox_genes <- c('Gzma','Gzmb','Gzmm','Gzmk','Prf1','Ctsw')
# seu_ys_nk <- AddModuleScore(seu_ys_nk, features = list(cytotox_genes), ctrl=5, name='cytotox_mod')
g <- VlnPlot(seu_ys_nk, features = 'cytotox_mod1', group.by = 'temp', idents = c('wt19_early','wt19_late','ko13_early','ko13_late','ko18_early','ko18_late'), cols=alpha(c('#8C4484','#8C4484','#ff8989','#ff8989','#DC0000FF','#DC0000FF'),0.8), adjust = 1.25) + 
  theme(plot.title = element_blank(), legend.position = 'none', axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank()) 
g$data$ident <- factor(g$data$ident, levels=c('wt19_early','wt19_late','ko13_early','ko13_late','ko18_early','ko18_late'))
print(g)

wilcox.test(seu_ys_nk@meta.data[which(seu_ys_nk@meta.data$temp == 'wt19_early'),'cytotox_mod1'], 
            seu_ys_nk@meta.data[which(seu_ys_nk@meta.data$temp == 'wt19_late'),'cytotox_mod1'])
wilcox.test(seu_ys_nk@meta.data[which(seu_ys_nk@meta.data$temp == 'ko13_early'),'cytotox_mod1'], 
            seu_ys_nk@meta.data[which(seu_ys_nk@meta.data$temp == 'ko13_late'),'cytotox_mod1'])
wilcox.test(seu_ys_nk@meta.data[which(seu_ys_nk@meta.data$temp == 'ko18_early'),'cytotox_mod1'], 
            seu_ys_nk@meta.data[which(seu_ys_nk@meta.data$temp == 'ko18_late'),'cytotox_mod1'])
