
#Figure4B
library(tidyverse)
library(ggsci)
gene <- readRDS("H:/Project/Smart/MMP-1/Figure2/lasso_fea.Rds")
GSE47460_anno_risk <- readRDS("H:/Project/Smart/MMP-1/Figure2/GSE47460_anno_risk.Rds")
imm_prop <- read.csv('estimation_matrix.csv',header = T,row.names = 1)
load("H:/Project/Smart/newmodel/runtime/Step0.database/GSE47460/GSE47460_expr.Rdata")
mycol <- pal_lancet()(2)

cibersort <- imm_prop[grep('SORT$',rownames(imm_prop),value = T),]
rownames(cibersort) <- str_split(rownames(cibersort), '_', simplify = T)[,1]

library(ggstatsplot)
library(ggpubr)
library(reshape2)
library(WGCNA)
colnames(imm_prop) == rownames(metadata)
Sample <- rownames(metadata[metadata$final.diagnosis == 'ILD',])

cibersort <- as.data.frame(t(cibersort[Sample]))
dat <- as.data.frame(t(dat[Sample]))[gene]
rownames(dat) == rownames(GSE47460_anno_risk)
dat$RiskScore = GSE47460_anno_risk$RiskScore
rownames(dat) == rownames(cibersort)

cor_pval <- corAndPvalue(cibersort,dat)
plot_dat <- melt(cor_pval$cor)

for (i in 1:nrow(plot_dat)) {
  plot_dat$Pvalue[i] <- cor_pval$p[plot_dat$Var1[i],plot_dat$Var2[i]]
}

p1 <- ggplot(plot_dat,aes(x=Var2,y=Var1)) #热图绘制
p2 <- p1+scale_fill_gradient2(low = mycol[1], mid = 'white', high = mycol[2],)+
  theme_minimal()+
  geom_point(aes(size=-log10(Pvalue),
                 fill=value), shape= 21)+
  theme(axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)
p2
ggsave(plot = p2,'cibersort_cor_IPF.pdf', width = 7, height = 7)

#Figure4C
rownames(dat) == rownames(cibersort)
corplot_dat <- cbind(cibersort, dat)
for (i in colnames(dat)) {
  for (j in colnames(cibersort)) {
    tmp <- corplot_dat[c(i,j)]
    colnames(tmp) <- c('X','Y')
    p <-  ggscatter(tmp,x = 'X', y = "Y",
                    color = "black", size = 3, # 点的颜色，大小
                    add = "reg.line",  # 添加回归线
                    add.params = list(color = "blue", fill = "lightgray"), # 回归线的调整
                    conf.int = TRUE, # 回归线的置信区间
                    cor.coef = TRUE, # 添加相关系数
                    cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                    ylab = j,
                    xlab = i)
    ggplot2::ggsave(plot = p, paste0(i,'_',j,'_IPF.pdf'), width = 5, height = 5)
  }
}

#Figure4D
genelist <- read.table('GeneList.txt', header = T,sep = '\t')
pathway <- split(genelist,genelist$Category)
pathway <- lapply(pathway, function(x) x <- x$Symbol)

library(clusterProfiler)
library(GSEABase)
library(GSVA)
library(ggpubr)
load("H:/Project/Smart/newmodel/runtime/Step0.database/GSE47460/GSE47460_expr.Rdata")
head(metadata)
metadata$Sex <- as.factor(metadata$Sex)
metadata$gold.stage <- as.factor(metadata$gold.stage)
metadata$smoker. <- as.factor(metadata$smoker.)
metadata$ild.subtype <- as.factor(metadata$ild.subtype)
metadata <- metadata[c(1,2,4:10,12,13)]
colnames(metadata) <- c('Age','Sex','Predicted.DLCO',"Predicted.FVC.pre.bd",
                        "Predicted.FVC.post.bd","Predicted.FEV1.pre.bd",
                        "Predicted.FEV1.post.bd","Emphysema.f950","Gold.stage",               
                        "Smoker","Diagnosis")
metadata$Diagnosis <- factor(metadata$Diagnosis, levels = c('Control','2-UIP/IPF'),
                             labels = c('Control','IPF'))

# dat <- expr_list$GSE47460
dat <- as.matrix(dat[Sample])
ssgsea_result <- gsva(dat, pathway, method = 'ssgsea')

metadata <- metadata[Sample,]
rownames(metadata) == rownames(GSE47460_anno_risk)
metadata$RiskScore = GSE47460_anno_risk$RiskScore


n <- t(scale(t(ssgsea_result[,Sample])))
n[n > 2] <- 2
n[n < -2] <- -2
rownames(metadata) == colnames(n)
pheatmap::pheatmap(n[,order(metadata$RiskScore)],annotation_col = metadata,
                   cluster_cols = T, show_colnames = F,
                   filename = 'pathway_imm_ipf.pdf',
                   width = 10, height = 5.5, cluster_rows = F)

#Figure4E
gene <- readRDS("H:/Project/Smart/MMP-1/Figure2/lasso_fea.Rds")
dat <- as.data.frame(t(dat))[gene]
rownames(dat) == rownames(GSE47460_anno_risk)
dat$RiskScore <- GSE47460_anno_risk$RiskScore
xcell <- t(ssgsea_result)
rownames(dat) == rownames(xcell)

cor_pval <- corAndPvalue(xcell,dat)
plot_dat <- melt(cor_pval$cor)

for (i in 1:nrow(plot_dat)) {
  plot_dat$Pvalue[i] <- cor_pval$p[plot_dat$Var1[i],plot_dat$Var2[i]]
}

p1 <- ggplot(plot_dat,aes(x=Var2,y=Var1)) #热图绘制
p2 <- p1+scale_fill_gradient2(low = mycol[1], mid = 'white', high = mycol[2],)+
  theme_minimal()+
  geom_point(aes(size=-log10(Pvalue),
                 fill=value), shape= 21)+
  theme(axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)
p2
ggsave(plot = p2,'mmp_cor_ipf.pdf', width = 8, height = 10)


#Figure4F

rownames(dat) == rownames(xcell)
corplot_dat <- cbind(xcell, dat)
for (i in colnames(dat)) {
  for (j in colnames(xcell)) {
    tmp <- corplot_dat[c(i,j)]
    colnames(tmp) <- c('X','Y')
    p <-  ggscatter(tmp,x = 'X', y = "Y",
                    color = "black", size = 3, # 点的颜色，大小
                    add = "reg.line",  # 添加回归线
                    add.params = list(color = "blue", fill = "lightgray"), # 回归线的调整
                    conf.int = TRUE, # 回归线的置信区间
                    cor.coef = TRUE, # 添加相关系数
                    cor.coeff.args = list(method = "pearson", label.sep = "\n"),
                    ylab = j,
                    xlab = i)
    ggplot2::ggsave(plot = p, paste0(i,'_',j,'_ipf.pdf'), width = 5, height = 5)
  }
}
