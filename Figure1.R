#Figure1A
#gene circle

library(RCircos)

# 绘制人染色体圈图 ------------------------------------------------------
# 导入内建人类染色体数据
data(UCSC.HG38.Human.CytoBandIdeogram)

# 设置染色体数据
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram
# 设置不显示的染色体，如 c(1,3)  
chr.exclude <- NULL
# 设置内部环形个数
tracks.inside <- 3
# 设置外部环形个数  
tracks.outside <- 0

# 导入上面四个基本参数
RCircos.Set.Core.Components(cyto.info, chr.exclude,tracks.inside, tracks.outside)

# 列出所有绘图参数
RCircos.List.Plot.Parameters()

# 绘制染色体图形，默认方法显示染色体名称。
RCircos.Set.Plot.Area()     
RCircos.Chromosome.Ideogram.Plot() 

# 添加基因名称与连线 ------------------------------------------------------

# 加载内置的RCircos.Gene.Label.Data数据集
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb=TxDb.Hsapiens.UCSC.hg38.knownGene
gene=genes(txdb)
gene <- as.data.frame(gene)

library(org.Hs.eg.db)
library(clusterProfiler)
df <- bitr(gene$gene_id, fromType = "ENTREZID",
           toType = c("SYMBOL"),
           OrgDb = org.Hs.eg.db)
head(df)
gene=merge(gene,df,by.y='ENTREZID',by.x='gene_id')

dat <- gene[gene$SYMBOL %in% grep('^MMP',gene$SYMBOL, value = T),c(2,3,4,7)]
colnames(dat) <- c("Chromosome","chromStart","chromEnd","Gene")


out.file <- "RCircosDemoHumanGenome1.pdf"
pdf(file=out.file, height=8, width=8, compress=TRUE)
RCircos.Set.Plot.Area()

par(mai=c(0.25, 0.25, 0.25, 0.25))
plot.new()
plot.window(c(-3,3), c(-3, 3))

RCircos.Chromosome.Ideogram.Plot()

side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(dat, track.num, side)

# 在染色体上添加基因名称， 指定内容在第几个环形生成
name.col <- 5
track.num <- 2
RCircos.Gene.Name.Plot(dat, name.col,track.num, side)

dev.off()


##Figure2B

library(tidyverse)
library(UpSetR)
library(ggsci)
library(scales)
mypal <- pal_lancet()(9)
mypal

show_col(mypal)


all_DEG <- readRDS("../data/all_DEG.Rds")
deg <- readRDS("../data/mmp_deg.Rds")

dat <- deg[grep('log',colnames(deg),value = T)]
colnames(dat) <- substr(colnames(dat),1,nchar(colnames(dat))-6)
dat$gene <- rownames(dat)

library(reshape2)
dat <- melt(dat,'gene')
dat[is.na(dat)] <- 0
for (i in 1:nrow(dat)) {
  dat$Pvalue.adjust[i] <- deg[dat$gene[i],paste0(dat$variable[i],'_adj')]
}
colnames(dat) <- c('Gene','Database','LogFC','Pvalue.adjust')

library(ggpubr)
p <- ggbarplot(dat,x = 'Gene',y = 'LogFC', fill = 'Database',palette = mypal[3:6],
               ylab = 'Log2FoldChange')+
  theme(axis.text.x = element_text(angle=60, hjust=1))+
  geom_hline(aes(yintercept = 0))


dat$LogFC[dat$LogFC < 0] <- -1
dat$LogFC[dat$LogFC > 0] <- 1
p1 <- ggplot(dat,aes(x=Gene,y=Database)) #热图绘制
p2 <- p1+scale_fill_gradient2(low = '#00468BFF', mid = 'white', high = '#ED0000FF',)+
  theme_bw()+
  geom_point(aes(size=-log10(Pvalue.adjust),
                 fill=LogFC), shape= 21)+
  theme(panel.grid = element_blank(),axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))
p2

library(patchwork)
p/p2


#Figure1C
ipf_expr <- readRDS("../data/ipf_expr.Rds")
library(corrplot)
for (i in names(ipf_expr)) {
  M <- cor(t(ipf_expr[[i]]))
  pdf(paste0(i,'_corplot.pdf'), width = 15, height = 15)
  corrplot.mixed(M)
  recordPlot()
  dev.off()
}

#Figure1E
mmp_expr <- readRDS("../data/mmp_expr.Rds")
dat <- as.data.frame(scale(t(mmp_expr[["GSE124685"]])))
library(tidyverse)
dat$Sample <- as.factor(str_split(rownames(dat),'_',simplify = T)[,1])
dat <- dat[sort(colnames(dat))]
library(reshape2)
dat <- melt(dat,'Sample')
library(ggpubr)  
ggboxplot(dat, x = 'Sample', y = 'value',
          palette = 'lancet', fill = 'Sample',
          facet.by = 'variable',short.panel.labs = T,nrow = 2,ggtheme = theme_bw())+
  stat_compare_means(label = 'p.signif')

ggsave('IPF1-3.pdf',width = 15, height = 6)