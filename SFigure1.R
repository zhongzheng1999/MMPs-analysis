
#SFigure1A
library(ggsci)
library(tidyverse)
library(FactoMineR)
load("H:/Project/Small/MMP/data/GSE32537/GSE32537_expr.Rdata")
gene <- sort(rownames(genes_expr)[grep('^MMP', rownames(genes_expr))])
dat <- as.data.frame(t(na.omit(genes_expr[gene,])))
metadata <- metadata[rownames(dat),]
gene.pca <- PCA(dat, ncp = 2, scale.unit = TRUE, graph = FALSE)

dat$State <- metadata$`final diagnosis:ch1`
#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)
pca_sample <- cbind(pca_sample, metadata[rownames(dat),])
pca_sample$samples <- rownames(pca_sample)
head(pca_sample)  #作图数据中包含了样本坐标和分组信息

pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = Group), size = 3) +  #根据样本坐标绘制二维散点图
  scale_color_lancet() +  #自定义颜色
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PC1:', pca_eig1, '%'), y = paste('PC2:', pca_eig2, '%'), color = '')+  #将 PCA 轴贡献度添加到坐标轴标题中
  stat_ellipse(aes(fill = Group), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  scale_fill_lancet()+
  theme_classic2()
p 
ggsave('GSE32537PCA.pdf',p,width = 6, height = 5)


library(reshape2)
dat <- melt(dat,'State')
colnames(dat)

library(ggpubr)
p1 <- ggboxplot(dat, x = "variable", y = "value",fill = "State", group.by = 'State',
                palette = 'lancet', ylab = 'Expression values', xlab = 'Matrix metalloproteinases')+
  stat_compare_means(aes(group = State),label = 'p.signif',hide.ns = T)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(plot = p1,filename = 'GSE32537_MMP_box.pdf',width = 9,height = 4.5)


#47460

load("H:/Project/Small/MMP/data/GSE47460/GSE47460_expr.Rdata")
gene <- sort(rownames(dat)[grep('^MMP', rownames(dat))])
dat <- as.data.frame(t(dat[gene,]))

metadata <- metadata[rownames(dat),]
gene.pca <- PCA(dat, ncp = 2, scale.unit = TRUE, graph = FALSE)

dat$State <- metadata$final.diagnosis
dat <- melt(dat,'State')
colnames(dat)

#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)
rownames(pca_sample) == rownames(metadata)
pca_sample <- cbind(pca_sample, metadata[rownames(pca_sample),])
pca_sample$samples <- rownames(pca_sample)
head(pca_sample)  #作图数据中包含了样本坐标和分组信息

pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = final.diagnosis), size = 3) +  #根据样本坐标绘制二维散点图
  scale_color_lancet() +  #自定义颜色
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PC1:', pca_eig1, '%'), y = paste('PC2:', pca_eig2, '%'), color = '')+  #将 PCA 轴贡献度添加到坐标轴标题中
  stat_ellipse(aes(fill = final.diagnosis), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  scale_fill_lancet()+
  theme_classic2()
p 
ggsave('GSE47460PCA.pdf',p,width = 6, height = 5)

p1 <- ggboxplot(dat, x = "variable", y = "value",fill = "State", group.by = 'State',
                palette = 'lancet', ylab = 'Expression values', xlab = 'Matrix metalloproteinases')+
  stat_compare_means(aes(group = State),label = 'p.signif',hide.ns = T)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(plot = p1,filename = 'GSE47460_MMP_box.pdf',width = 9,height = 4.5)


#GSE124685
id2name <- readRDS("H:/Project/main/main/gtf/hsv27_id2name.rds")
fpkms <- readRDS("H:/Project/Small/MMP/data/GSE124685/GSE124685fpkms.Rds")
mmp_index <- sort(grep('^MMP',id2name$gene_name,value = T))

id2name <- id2name[id2name$gene_name %in% mmp_index,]
dat <- cbind(fpkms[id2name$gene_id,],id2name$gene_name)
rownames(dat) <- dat$`id2name$gene_name`
dat <- dat[-ncol(dat)]

dat <- as.data.frame(t(na.omit(dat)))

##pca
gene.pca <- PCA(dat, ncp = 2, scale.unit = TRUE, graph = FALSE)

#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)
pca_sample$samples <- ifelse(rownames(pca_sample) %in% grep('^CON',rownames(pca_sample),value = T),'Control','IPF')
head(pca_sample)  #作图数据中包含了样本坐标和分组信息

pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2],2 )

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = samples), size = 3) +  #根据样本坐标绘制二维散点图
  scale_color_lancet() +  #自定义颜色
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PC1:', pca_eig1, '%'), y = paste('PC2:', pca_eig2, '%'), color = '')+  #将 PCA 轴贡献度添加到坐标轴标题中
  stat_ellipse(aes(fill = samples), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  scale_fill_lancet()+
  theme_classic2()
p 
ggsave('GSE124685PCA.pdf',p,width = 6, height = 5)
###END

dat <- dat[sort(colnames(dat))]
dat <- as.data.frame(scale(dat))
dat$State <- ifelse(rownames(dat) %in% grep('^CON',rownames(dat),value = T),'Control','IPF')
dat <- melt(dat,'State')
colnames(dat)

p1 <- ggboxplot(dat, x = "variable", y = "value",fill = "State", group.by = 'State',
                palette = 'lancet', ylab = 'Scaled FPKM values', xlab = 'Matrix metalloproteinases')+
  stat_compare_means(aes(group = State),label = 'p.signif',hide.ns = T)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p1

ggsave(plot = p1,filename = 'GSE124685_MMP_box.pdf',width = 9,height = 4.5)


#GSE150910
id2name <- readRDS("H:/Project/main/main/gtf/hsv27_id2name.rds")
fpkms <- readRDS("H:/Project/Small/MMP/data/GSE150910/GSE150910fpkms.Rds")
mmp_index <- sort(grep('^MMP',id2name$gene_name,value = T))

id2name <- id2name[id2name$gene_name %in% mmp_index,]
dat <- cbind(fpkms[id2name$gene_name,],id2name$gene_name)
rownames(dat) <- dat$`id2name$gene_name`
dat <- dat[-ncol(dat)]
dat <- as.data.frame(t(na.omit(dat)))
dat <- dat[sort(colnames(dat))]
dat <- dat[c(grep('^ipf',rownames(dat),value = T),grep('^con',rownames(dat),value = T)),]
dat <- dat
##pca
gene.pca <- PCA(dat, ncp = 2, scale.unit = TRUE, graph = FALSE)
#提取样本在 PCA 前两轴中的坐标
pca_sample <- data.frame(gene.pca$ind$coord[ ,1:2])
head(pca_sample)
pca_sample$samples <- ifelse(rownames(pca_sample) %in% grep('^con',rownames(pca_sample),value = T),'Control','IPF')
head(pca_sample)  #作图数据中包含了样本坐标和分组信息

pca_eig1 <- round(gene.pca$eig[1,2], 2)
pca_eig2 <- round(gene.pca$eig[2,2], 2)

p <- ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = samples), size = 3) +  #根据样本坐标绘制二维散点图
  scale_color_lancet() +  #自定义颜色
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) +  #去除背景和网格线
  labs(x =  paste('PC1:', pca_eig1, '%'), y = paste('PC2:', pca_eig2, '%'), color = '')+  #将 PCA 轴贡献度添加到坐标轴标题中
  stat_ellipse(aes(fill = samples), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +
  scale_fill_lancet()+
  theme_classic2()
p 
ggsave('GSE150910PCA.pdf',p,width = 6, height = 5)
###END

#SFigure1B
ipf_expr <- readRDS("../data/ipf_expr.Rds")
library(corrplot)
for (i in names(ipf_expr)) {
  M <- cor(t(ipf_expr[[i]]))
  pdf(paste0(i,'_corplot.pdf'), width = 15, height = 15)
  corrplot.mixed(M)
  recordPlot()
  dev.off()
}


#SFigure1C
library(WGCNA)
ipf_cor <- lapply(ipf_expr, function(x){
  x <- t(x)
  x <- corAndPvalue(x, x[,'MMP7'])
  x <- as.data.frame(cbind(x[[1]],x[[2]]))
  colnames(x) <- c('Cor','Pvalue')
  x$Pvalue.adj <- p.adjust(x$Pvalue,method = 'BH')
  x$gene <- rownames(x)
  x$pair <- paste0('MMP7-',x$gene)
  return(x)
})

up_pathway <- lapply(ipf_cor, function(x){
  x <- x[x$Cor > 0 & x$Pvalue < 0.05,]
  return(x$pair)
})
down_pathway <- lapply(ipf_cor, function(x){
  x <- x[x$Cor < 0 & x$Pvalue < 0.05,]
  return(x$pair)
})

library(VennDiagram)
venn.plot <- venn.diagram(up_pathway,col="white",fill = ggsci::pal_lancet()(4),
                          filename=NULL,main = 'Positive Correlation')

#将venn.plot通过grid.draw画到pdf文件中
pdf("pathway_up_venn.pdf", width = 4, height = 4)
grid.draw(venn.plot)
file.remove(list.files(pattern = '*log'))
dev.off()


venn.plot <- venn.diagram(down_pathway,col="white",fill = ggsci::pal_lancet()(4),
                          filename=NULL,main = 'Negative Correlation')

#将venn.plot通过grid.draw画到pdf文件中
pdf("pathway_down_venn.pdf", width = 4, height = 4)
grid.draw(venn.plot)
file.remove(list.files(pattern = '*log'))
dev.off()

up_pathway_ident <- Reduce(intersect,up_pathway)
down_pathway_ident <- Reduce(intersect,down_pathway)

#SFigure1D
GSE150910fpkms <- readRDS("../data/GSE150910fpkms.Rds")
dat <- GSE150910fpkms[sort(grep('^MMP',rownames(GSE150910fpkms),value = T)),]
dat <- as.data.frame(scale(t(dat)))

library(tidyverse)
dat$group <- str_split(rownames(dat), '_', simplify = T)[,1]
dat$group

library(reshape2)
dat_v <- melt(dat,c('group'))
head(dat_v)
dat_v$group <- factor(dat_v$group, levels = c('control','chp','ipf'), labels = c('Control','CHP','IPF'))


library(ggpubr)
library(ggsci)
library(scales)
mycol <- pal_lancet()(9)
show_col(mycol)

ggboxplot(dat_v, x = 'variable', y = 'value', fill = 'group', palette = mycol[c(1,3,2)])+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6))
ggsave('GSE150910_CHP.pdf', width = 11, height = 5)

dat <- t(GSE150910fpkms[sort(grep('^MMP',rownames(GSE150910fpkms),value = T)),])
pca <- prcomp(dat, scale. = T)
pca_dat <- data.frame(MMP_Score = pca$x[,1])
pca_dat$group <- str_split(rownames(dat), '_', simplify = T)[,1]
pca_dat$group <- factor(pca_dat$group, levels = c('control','chp','ipf'), labels = c('Control','CHP','IPF'))

p3 <- pca_dat %>% ggplot(aes(x = MMP_Score)) +
  stat_ecdf(aes(color = group),
            geom = "smooth",
            size = 1.5) +
  scale_color_manual(values = mycol[c(1,3,2)]) +
  labs(y = "ECDF",title = 'MMP Family Empirical Cumulative Distribution')+
  theme_bw()
p3

ggsave('MMP_ECDF.pdf',height = 5,width = 6)


#SFigure1E
dat <- openxlsx::read.xlsx('../data/GSE124692.xlsx', sheet = 1, rowNames = T)
meta <- openxlsx::read.xlsx('../data/GSE124692.xlsx', sheet = 2, rowNames = T)
geneid <- openxlsx::read.xlsx('../data/GSE124692.xlsx', sheet = 3, rowNames = T)

rownames(meta) == colnames(dat)
colnames(dat) <- meta$SampleName

geneid <- geneid[rownames(dat),]
geneid <- geneid[!duplicated(geneid$GeneName),]
dat <- dat[rownames(geneid),]
table(rownames(dat) == rownames(geneid))
rownames(dat) <- geneid$GeneName

dat_mmp <- as.data.frame(t(dat[sort(grep('^MMP',rownames(dat),value = T)),]))
# dat <- as.data.frame(scale(t(dat)))

library(tidyverse)
dat_mmp$group <- str_split(rownames(dat_mmp), '_', simplify = T)[,1]
dat_mmp$group

library(reshape2)

dat_v <- melt(dat_mmp,c('group'))
head(dat_v)
dat_v$group <- ifelse(dat_v$group == 'normal', 'Normal', dat_v$group)
dat_v$group <- factor(dat_v$group, levels = c('Normal','ALI','IPF'), labels = c('Normal','ALI','IPF'))


library(ggpubr)
library(ggsci)
library(scales)
mycol <- pal_lancet()(9)
show_col(mycol)

ggboxplot(dat_v, x = 'variable', y = 'value', fill = 'group', palette = mycol[c(1,3,2)])+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.6))
ggsave('GSE134692_ALI.pdf', width = 11, height = 5)

dat_mmp <- t(dat[sort(grep('^MMP',rownames(dat),value = T)),])
pca <- prcomp(dat_mmp, scale. = T)
pca_dat <- data.frame(MMP_Score = pca$x[,1])
pca_dat$group <- str_split(rownames(dat_mmp), '_', simplify = T)[,1]
pca_dat$group <- ifelse(pca_dat$group == 'normal', 'Normal', pca_dat$group)
pca_dat$group <- factor(pca_dat$group, levels = c('Normal','ALI','IPF'), labels = c('Normal','ALI','IPF'))

p3 <- pca_dat %>% ggplot(aes(x = MMP_Score)) +
  stat_ecdf(aes(color = group),
            geom = "smooth",
            size = 1.5) +
  scale_color_manual(values = mycol[c(1,3,2)]) +
  labs(y = "ECDF",title = 'MMP Family Empirical Cumulative Distribution')+
  theme_bw()
p3

ggsave('MMP_ECDF_GSE134692.pdf',height = 5, width = 6)

#SFigure1F
library(ggpubr)
mmp_expr <- readRDS("H:/Project/Smart/MMP/mmp/mmp_expr.Rds")
meta <- readRDS("H:/Project/Smart/MMP/mmp/GSE150910meta.Rds")
meta <- meta[meta$diagnosis == 'ipf',]

dat <- t(mmp_expr$GSE150910)
dat <- as.data.frame(scale(dat[rownames(meta),]))
rownames(dat) == rownames(meta)

dat$State <- meta$diagnosis
dat$rs35705950_genotype <- meta$rs35705950_genotype

library(reshape2)
dat <- melt(dat, c('State','rs35705950_genotype'))
dat$rs35705950_genotype <- as.factor(dat$rs35705950_genotype)
dat_log <- dat

ggboxplot(dat_log,x = 'rs35705950_genotype', y = 'value', facet.by = 'variable',
          fill = 'rs35705950_genotype', palette = 'lancet',nrow = 2,
          ylab = 'logNormalization Value', ggtheme = theme_bw())+
  stat_compare_means(label = 'p.signif')
ggsave('ipf_eQTL.pdf', height = 6, width = 16)