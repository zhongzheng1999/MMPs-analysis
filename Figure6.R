library(tidyverse)
library(Seurat)

#Figure6A
load("H:/Project/Smart/MMP/mmp/SingleCell/sce.Rdata")
DimPlot(sce_sub, split.by = 'Diagnosis')

macro <- subset(sce_sub, subset = celltype == "Macrophages")

gene <- rownames(macro)
mmp_gene <- grep('^MMP', gene, value = T)
mmp_gene

macro <- FindVariableFeatures(macro)
grep('^MMP', VariableFeatures(macro), value = T)
macro <- NormalizeData(macro)
macro <- ScaleData(macro)

macro <- SCTransform(macro)
macro <- RunPCA(macro)
macro <- RunUMAP(macro, dims = 1:30)
macro <- FindNeighbors(macro, dims = 1:30)
macro <- FindClusters(macro)
DimPlot(macro, label = TRUE) + NoLegend()
DotPlot(macro, features = grep('^MMP', VariableFeatures(macro), value = T))
FeaturePlot(macro, features = grep('^MMP', VariableFeatures(macro), value = T), max.cutoff = 1)


#high and low
DefaultAssay(macro) <- 'RNA'
highCells=colnames(subset(x = macro, subset = MMP19 > 0, slot = 'counts'))
highORlow=ifelse(colnames(macro) %in% highCells,'high','low')
table(highORlow)
macro@meta.data$highORlow=highORlow
markers <- FindMarkers(macro, ident.1 = "high", 
                       group.by = 'highORlow', 
                       subset.ident = "0",logfc.threshold = 0)
head(x = markers)

deg <- markers
deg[which(deg$p_val_adj < 0.01 & deg$avg_log2FC <= -2),'sig'] <- 'Down'
deg[which(deg$p_val_adj < 0.01 & deg$avg_log2FC>= 2),'sig'] <- 'Up'
deg[which(deg$p_val_adj >= 0.01 | abs(deg$avg_log2FC) < 2),'sig'] <- 'None'
sub = paste0('Up:',length(which((deg[,6]=='Up'))),' ', 'Down:',length(which((deg[,6]=='Down'))))

ggplot(deg, aes(x = avg_log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.6, size = 1) + #alpha调节透明度，size是点的大小
  scale_colour_manual(values  = c('red2', 'blue2', 'gray'), limits = c('Up', 'Down', 'None')) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), plot.title = element_text(hjust = 0.5)) +
  theme(legend.key = element_rect(fill = 'transparent'), legend.background = element_rect(fill = 'transparent'), legend.position = c(0.9, 0.93)) +
  geom_vline(xintercept = c(-2, 2), color = 'gray', size = 0.3) +
  geom_hline(yintercept = -log(0.01, 10), color = 'gray', size = 0.3) +
  #
  labs(x = '\nLog2 Fold Change', y = '-Log10(adj.P.Val)\n', color = '', title = 'MMP19+ Vs MMP19- (Macrophages)',subtitle = sub)+
  xlim(c(-8,8))
ggsave('Macrophages_volcano.pdf', width = 5, height = 6)

#Figure6B
##GOKEGG
#富集过程
library(clusterProfiler)
library(org.Hs.eg.db)
deg$gene_name <- rownames(deg)
df <- bitr(unique(deg$gene_name), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
deg=merge(deg,df,by.y='SYMBOL',by.x='gene_name')
gene_diff= deg[deg$sig != 'None','ENTREZID'] 
gene_all=as.character(deg[ ,'ENTREZID'] )

R.utils::setOption( "clusterProfiler.download.method",'auto' )
kk.result <- enrichKEGG(gene         = gene_diff,
                        organism     = 'hsa',
                        universe     = gene_all,
                        pvalueCutoff = 0.9,
                        qvalueCutoff =0.9,)
kk.result <- setReadable(kk.result, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

kk.result@result$'Enrichment Fold'=apply(kk.result@result,1,function(x){
  GeneRatio=eval(parse(text=x["GeneRatio"]))
  BgRatio=eval(parse(text=x["BgRatio"]))
  enrichment_fold=round(GeneRatio/BgRatio,2)
  return(enrichment_fold)
})


go_enrich_results <- lapply( c('BP','MF','CC') , function(ont) {
  ego <- enrichGO(gene          = gene_diff,
                  universe      = gene_all,
                  OrgDb         = org.Hs.eg.db,
                  ont           = ont ,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.99,
                  qvalueCutoff  = 0.99,
                  readable      = TRUE)
  
  ego@result$'Enrichment Fold'=apply(ego@result,1,function(x){
    GeneRatio=eval(parse(text=x["GeneRatio"]))
    BgRatio=eval(parse(text=x["BgRatio"]))
    enrichment_fold=round(GeneRatio/BgRatio,2)
    return(enrichment_fold)
  })
  #print( head(ego) )
  return(ego)
})


#富集结束，文本结果整理
dat_enrich <- list(GOBP = go_enrich_results[[1]],
                   GOMF = go_enrich_results[[2]],
                   GOCC = go_enrich_results[[3]],
                   KEGG = kk.result)

#文本结果下载
write.xlsx(dat_enrich,'enrich_results.xlsx')


#统一作图
enrich_plot_enrich <- list()
for ( i in names(dat_enrich)) {
  dat <- dat_enrich[[i]]@result %>% slice_min(pvalue,n=8)
  
  enrich_plot_enrich[[i]] <- ggplot(dat,aes(`Enrichment Fold`, fct_reorder(Description, `Enrichment Fold`)))+
    geom_segment(aes(xend=0, yend = Description))+
    geom_point(aes(color=p.adjust, size = Count)) +
    scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
    scale_size_continuous(range=c(2, 10)) +
    theme_classic(base_size = 16) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 40))+
    ylab(NULL)+ggtitle(paste(i))
}

#展示上调基因富集图
p3 <- enrich_plot_enrich[[1]]+enrich_plot_enrich[[2]]+enrich_plot_enrich[[3]]+enrich_plot_enrich[[4]]
ggsave(plot = p3,filename = 'Enrich_plot.pdf',width = 18,height = 10)

#Figure6C
library(CellChat)
library(tidyverse)
library(ggalluvial)
library(Seurat)

load("H:/Project/Smart/MMP/mmp/SingleCell/sce.Rdata")
DefaultAssay(sce_sub) <- 'RNA'
sce_sub <- NormalizeData(sce_sub)

macro <- subset(sce_sub, subset = celltype == "Macrophages")
DefaultAssay(macro) <- 'RNA'
highCells=colnames(subset(x = macro, subset = MMP19 > 0, slot = 'counts'))
highORlow=ifelse(colnames(macro) %in% highCells,'high','low')
table(highORlow)
macro@meta.data$highORlow=highORlow
macro@meta.data[["highORlow"]] <- paste0('MMP19 ',macro@meta.data[["highORlow"]])

meta <- data.frame(barcode = rownames(macro@meta.data), celltype = macro$highORlow)

sce_sub@meta.data[meta$barcode,"celltype"] <- meta$celltype


cellchat <- createCellChat(object = sce_sub,group.by = "celltype")

CellChatDB <- CellChatDB.human
# CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, population.size = T)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
# cellchat <- computeNetSimilarity(cellchat, type = "functional")
# cellchat <- netEmbedding(cellchat, type = "functional")
# cellchat <- netClustering(cellchat, type = "functional")

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight

mmp_weight <- as.data.frame(t(mat[c('MMP19 high', 'MMP19 low'),]))
mmp_weight$high.scale <- (mmp_weight$`MMP19 high`-min(mmp_weight$`MMP19 high`))/(max(mmp_weight$`MMP19 high`)-min(mmp_weight$`MMP19 high`))
mmp_weight$low.scale <- (mmp_weight$`MMP19 low`-min(mmp_weight$`MMP19 low`))/(max(mmp_weight$`MMP19 low`)-min(mmp_weight$`MMP19 low`))
mmp_weight$CellType <- rownames(mmp_weight)
mmp_weight <- mmp_weight[order(mmp_weight$high.scale),]
mmp_weight$CellType <- factor(mmp_weight$CellType,levels = mmp_weight$CellType)
p1 <- ggplot(mmp_weight, aes(CellType, high.scale))+
  geom_bar(stat="identity", fill = "red")+
  ylab(c('Normalization Weight'))+xlab(NULL)+
  theme_bw()+
  theme(axis.text.x = element_blank())
p2 <- ggplot(mmp_weight, aes(CellType, low.scale))+
  geom_bar(stat="identity", fill = "blue")+
  ylab(c('Normalization Weight'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 35,vjust = 1,hjust = 1))
library(patchwork)
p1/p2

ggsave('weight.pdf', height = 6, width = 10)

#Figure6D
mmp_weight$rank.high <- rank(mmp_weight$high.scale)
mmp_weight$rank.low <- rank(mmp_weight$low.scale)

dat <- mmp_weight[5:7]
tbody.style = tbody_style(color = "black", hjust=1, x=0.9)
ggtexttable(dat, rows = NULL,
            theme = ttheme(
              colnames.style = colnames_style(color = "black", fill = "grey"),
              tbody.style = tbody.style
            ))


netVisual_bubble(cellchat, sources.use = c(15,16), targets.use = c(2,6,8), remove.isolate = FALSE)
ggsave('function.pdf', width = 4, height = 8.5)

#Figure6E
library(bigSCale)
library(SingleCellExperiment)
library(Seurat)
library(igraph)
library(tidyverse)
AT2 <- subset(sce_sub, subset = celltype == "AT2")
sub_AT2 <- SplitObject(AT2, split.by = "Diagnosis")

expr.ctl <- sub_AT2$Control@assays$RNA@counts
expr.t2d <- sub_AT2$IPF@assays$RNA@counts
model=compute.network.model(expr.data = cbind(expr.ctl,expr.t2d))

gene.names <- rownames(expr.ctl)

results.ctl=compute.network(expr.data = expr.ctl,gene.names = gene.names,model = model, quantile.p = 0.6)
results.t2d=compute.network(expr.data = expr.t2d,gene.names = gene.names,model = model, quantile.p = 0.6)


output=homogenize.networks(list(results.ctl,results.t2d))
results.ctl=output[[1]]
results.t2d=output[[2]]
saveRDS(output,'output.Rds')

comparison=compare.centrality(list(results.ctl$centrality,results.t2d$centrality),c('Control','IPF'))
DT::datatable(comparison$PAGErank)

cor.ctl <- as.data.frame(results.ctl$correlations)
TNFRSF12A.ctl <- cor.ctl[,'TNFRSF12A']
names(TNFRSF12A.ctl) <- rownames(cor.ctl)
TNFRSF12A.ctl[TNFRSF12A.ctl > 0.6]
TNFRSF12A.ctl[TNFRSF12A.ctl < -0.6]

cor.t2d <- as.data.frame(results.t2d$correlations)
TNFRSF12A.t2d <- cor.t2d[,'TNFRSF12A']
names(TNFRSF12A.t2d) <- rownames(cor.t2d)
TNFRSF12A.t2d[TNFRSF12A.t2d > 0.6]
TNFRSF12A.t2d[TNFRSF12A.t2d < -0.6]

gene <- names(TNFRSF12A.t2d[TNFRSF12A.t2d > 0.6])
dat <- data.frame(TNFRSF12A.t2d[TNFRSF12A.t2d > 0.6])
dat$S <- 'TNFRSF12A'
dat$R <- rownames(dat)
colnames(dat) <- c('Cor','S','R')
write.csv(dat, 'TNFRSF12A_IPF.csv')
# toCytoscape(results.ctl$graph,'results.ctl.json')
# toCytoscape(results.t2d$graph,'results.t2d.json')
