library(Seurat)

#SFigure6A
load("H:/Project/Smart/MMP/mmp/SingleCell/sce.Rdata")
DimPlot(sce_sub)
FeaturePlot(sce_sub, features = 'TNFRSF12A', max.cutoff = 1,cols = c('grey','red'))
ggsave( filename = 'TNFRSF12A.pdf',width = 6,height = 5)

#SFigure6B
#Refer to the "IPF cell atlas" web

#SFigure6C
library(corrplot)
cor.dat.ctl <- cor.ctl[intersect(gene, rownames(cor.ctl)),intersect(gene, colnames(cor.ctl))]
corrplot(as.matrix(cor.dat.ctl))

cor.dat.t2d <- cor.t2d[intersect(gene, rownames(cor.ctl)),intersect(gene, colnames(cor.ctl))]
corrplot(as.matrix(cor.dat.t2d))

dat <- cor.dat.ctl
dat[lower.tri(dat)] <- cor.dat.t2d[lower.tri(cor.dat.t2d)]

bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
pheatmap::pheatmap(dat, cluster_rows = F, cluster_cols = F, 
                   color = colorRampPalette(colors = c("#000080","white","red"))(100),
                   legend_breaks=seq(-1,1,0.3),
                   show_rownames = F, show_colnames = F)

#SFigure6D
##GOKEGG
#富集过程
library(clusterProfiler)
library(org.Hs.eg.db)
deg <- data.frame(gene_name = rownames(cor.t2d), row.names = rownames(cor.t2d))
deg$sig <- 'None'
deg[names(TNFRSF12A.t2d[TNFRSF12A.t2d > 0.6]),'sig'] <- 'pos'
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
openxlsx::write.xlsx(dat_enrich,'enrich_results.xlsx')


#统一作图
enrich_plot_enrich <- list()
for ( i in names(dat_enrich)) {
  dat <- dat_enrich[[i]]@result %>% slice_min(pvalue,n=10)
  
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
ggsave(plot = p3,filename = 'Enrich_plot.pdf',width = 18,height = 6)
