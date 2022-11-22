library(Seurat)
load("H:/Project/Smart/MMP/mmp/SingleCell/sce.Rdata")
DimPlot(sce_sub)
sce_sub <- FindVariableFeatures(sce_sub)
sce_sub <- NormalizeData(sce_sub)
sce_sub <- ScaleData(sce_sub)
marker <- FindAllMarkers(sce_sub, min.pct = 0.5, only.pos = T)
# saveRDS(marker,'marker.Rds')
marker_selected <- markers  %>% 
  dplyr::filter(p_val_adj < 0.05) %>% 
  dplyr::filter(pct.1 >= 0.7 & pct.2 <= 0.5) %>% 
  dplyr::group_by(cluster) %>% 
  dplyr::slice_max(order_by = avg_log2FC, n = 1)
DotPlot(sce_sub,features = marker_selected$gene)

#SFigure5B
cells <- FetchData(sce, 
                   vars = c("ident", "Sample")) %>%
  dplyr::count(ident, Sample) %>%
  tidyr::spread(ident, n)
cells[is.na(cells)] <- 0
cells_prop <- as.data.frame(prop.table(as.matrix(cells[2:15]),margin = 1))
cells_prop$Sample <- cells$Sample
write.xlsx(list(cells=cells,cells_prop = cells_prop),'All_Sample_cells.xlsx',col.names = T)

library(ggpubr)
library(reshape2)
dat <- melt(cells_prop,'Sample')
dat$State <- ifelse(substr(dat$Sample,1,1) == 'C','Healthy','Diseased' )
dat$State <- as.factor(dat$State) %>% relevel('Healthy')


p1 <- ggboxplot(dat, x = "variable", y = "value",fill = "State", group.by = 'State',
                palette = 'lancet', ylab = 'Cell Proportion', xlab = '')+
  stat_compare_means(aes(group = State),label = 'p.signif',hide.ns = T)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1
ggsave(plot = p1,filename = 'Celltype_State_box.pdf',width = 9,height = 6)

#Figure5C
#Refer to Figure3

#Figure5D
library(AUCell)
library(GSEABase)
gmts <- list.files('H:/Project/Smart/newmodel/Msigdb',pattern = 'gmt', full.names = T)
geneset <- getGmt(gmts[4])
dat <- sce_sub@assays$RNA@counts
geneset <- subsetGeneSets(geneSets = geneset, rownames(dat))
geneset <- setGeneSetNames(geneSets = geneset,newNames = paste(names(geneset),"(",nGenes(geneset),"g)",sep = ""))
set.seed(321)
cells_rankings <- AUCell_buildRankings(dat,plotStats = T,splitByBlocks=TRUE)  # 关键一步
cells_AUC <- AUCell_calcAUC(geneset, cells_rankings)#, aucMaxRank=nrow(cells_rankings)*0.08
save(cells_AUC,file = 'cells_AUC.Rdata')
ident_num <- which(str_split(cells_AUC@NAMES, "\\(", simplify = T)[,1] %in% ident_pathway)
cells_assignment <- AUCell_exploreThresholds(cells_AUC[1,], plotHist=TRUE,assign=TRUE,nBreaks = 150)#
# ggsave('PyroptosisAUC_Thresholds.pdf',height = 5,width = 7)
save(cells_assignment,file = 'cells_assignment.Rdata')
cells_assignment$geneSet$aucThr$thresholds
cells_assignment$geneSet$aucThr$selected
aucs <- t(getAUC(cells_AUC))
aucs_bi <- ifelse(aucs > cells_assignment$geneSet$aucThr$selected, 1,0)
colnames(aucs) <- str_split(colnames(aucs), "\\(", simplify = T)[,1]
table(colnames(sce_sub) == rownames(aucs))
sce_sub <- AddMetaData(sce_sub,metadata = aucs, col.name = colnames(aucs))

for (i in colnames(aucs)) {
  p <- FeaturePlot(sce_sub,features = i,reduction = 'umap', cols = c('grey','red'))
  ggsave(plot = p, filename = paste(i,'.png'),width = 6,height = 5)
}

FeaturePlot(sce_sub, features = c('MMP19', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB'), blend = T, blend.threshold = 0.1)
ggsave( filename = 'MMP19_pathway.pdf',width = 14,height = 5)
