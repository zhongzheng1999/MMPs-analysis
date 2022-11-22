library(glmnet)
fit <- readRDS("H:/Project/Smart/MMP-1/Figure2/result_fit.Rds")
lasso_fea <- readRDS("H:/Project/Smart/MMP-1/Figure2/lasso_fea.Rds")
GSE47460 <- readRDS("H:/Project/Smart/MMP-1/Figure2/GSE47460_anno_risk.Rds")
GSE32537 <- readRDS("H:/Project/Smart/MMP-1/Figure2/GSE32537_anno_risk.Rds")
ipf_expr <- readRDS("H:/Project/Smart/MMP-1/ipf_expr.Rds")


###Figure3AB
anno <- GSE47460#GSE32537

coef_dat <- data.frame(Pvalue = 0, Coef = 0, confint.lower = 0, confint.upper = 0)

age <- lm(anno$meta_lung_function ~ anno$age)
tmp <- summary(age)
tmp$coefficients[2,4]
coef_dat$Coef[1] <- signif(coef(age)[2],2)
coef_dat$confint.lower[1] <- signif(confint(age)[2,1],2)
coef_dat$confint.upper[1] <- signif(confint(age)[2,2],2)
coef_dat$Pvalue[1] <- signif(tmp$coefficients[2,4],2)

sex <- lm(anno$meta_lung_function ~ anno$gender)
tmp <- summary(sex)
coef_dat[2,] <- NA
coef_dat$Coef[2] <- signif(coef(sex)[2],2)
coef_dat$confint.lower[2] <- signif(confint(sex)[2,1],2)
coef_dat$confint.upper[2] <- signif(confint(sex)[2,2],2)
coef_dat$Pvalue[2] <- signif(tmp$coefficients[2,4],2)

func <- lm(anno$meta_lung_function ~ anno$RiskScore) 
tmp <- summary(func)
coef_dat[3,] <- NA
coef_dat$Coef[3] <- signif(coef(func)[2],2)
coef_dat$confint.lower[3] <- signif(confint(func)[2,1],2)
coef_dat$confint.upper[3] <- signif(confint(func)[2,2],2)
coef_dat$Pvalue[3] <- signif(tmp$coefficients[2,4],2)

library(forestplot)
mydata <- coef_dat
rownames(mydata) <-  c("Age","Gender","RiskScore") 
mydata$Coef95CI <- paste0(mydata$Coef,"(",mydata$confint.lower,' - ',mydata$confint.upper,')')
mydata$factor <- rownames(mydata)
mydata <- mydata[c(6,1,5,2:4)]

pdf('univariate_line_GSE47460.pdf',width = 10,height = 4,onefile = F)
forestplot(labeltext=as.matrix(mydata[,c(1:3)]),#只展示前面的三列
           mean= mydata$Coef,#OR值
           lower= mydata$confint.lower,#95%CI下限
           upper= mydata$confint.upper,#95%CI上限
           zero=0,#OR值的位置，如果是线性回归则写0
           boxsize=0.2,#中间方框的大小
           # xticks=seq(1,50,190),#x轴刻度
           lwd.zero=2,#中间竖线的宽度
           lwd.ci=2, 
           col=fpColors(box='black',lines = 'black',zero = 'gray'),#颜色，box，lines和zero分别是方框，线条，中间竖线的颜色
           xlab="Coeffecients",#x轴标签
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 1.5),xlab  = gpar(cex = 1.5),
                            cex = 1.5),#设置字体大小
           lty.ci = "solid",
           title = "Forestplot", #标题
           graph.pos=2#中间竖线的位置
)
dev.off()

###Figure3CD
all <- lm(anno$meta_lung_function ~ anno$age+anno$gender+anno$RiskScore) 
tmp <- summary(all)
tmp <- data.frame(tmp$coefficients,confint(all))
tmp <- tmp[-1,c(4,1,5,6)]

library(forestplot)
mydata <- tmp
mydata <- apply(mydata, 2, function(x){
  x <- signif(x,2)
})
rownames(mydata) <-  c("Age","Gender","RiskScore")
mydata <- as.data.frame(mydata)
mydata$Coef95CI <- paste0(mydata$Estimate,"(",mydata$X2.5..,' - ',mydata$X97.5..,')')
mydata$factor <- rownames(mydata)
mydata <- mydata[c(6,1,5,2:4)]

pdf('multivariate_line_GSE47460.pdf',width = 10,height = 4,onefile = F)
forestplot(labeltext=as.matrix(mydata[,c(1:3)]),#只展示前面的三列
           mean= mydata$Estimate,#OR值
           lower= mydata$X2.5..,#95%CI下限
           upper= mydata$X97.5..,#95%CI上限
           zero=0,#OR值的位置，如果是线性回归则写0
           boxsize=0.2,#中间方框的大小
           # xticks=seq(1,50,190),#x轴刻度
           lwd.zero=2,#中间竖线的宽度
           lwd.ci=2, 
           col=fpColors(box='black',lines = 'black',zero = 'gray'),#颜色，box，lines和zero分别是方框，线条，中间竖线的颜色
           xlab="Coeffecients",#x轴标签
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 1.5),xlab  = gpar(cex = 1.5),
                            cex = 1.5),#设置字体大小
           lty.ci = "solid",
           title = "Forestplot", #标题
           graph.pos=2#中间竖线的位置
)
dev.off()

#Figure3EF
lasso_fea <- readRDS("H:/Project/Smart/MMP-1/Figure2/lasso_fea.Rds")
GSE47460 <- readRDS("H:/Project/Smart/MMP-1/Figure2/GSE47460_anno_risk.Rds")
GSE32537 <- readRDS("H:/Project/Smart/MMP-1/Figure2/GSE32537_anno_risk.Rds")
ipf_expr <- readRDS("H:/Project/Smart/MMP-1/ipf_expr.Rds")
anno <- GSE47460
dat <- ipf_expr$GSE47460[lasso_fea,]
colnames(dat) == rownames(anno)
n <- t(scale(t(dat)))
n[n > 2] = 2
n[n < -2] = -2
n <- na.omit(n)
pheatmap::pheatmap(n[,order(anno$RiskScore)], annotation_col = anno,
                   show_colnames = F, cluster_cols = F,
                   filename = 'pheamap_modelmmp_GSE47460.pdf',
                   width = 8.5, height = 5)

#Figure3G
library(WGCNA)
expr_list <- readRDS("H:/Project/Smart/MMP/data/expr_list.Rds")
for (i in names(expr_list)) {
  anno <- readRDS(paste0("H:/Project/Smart/MMP-1/Figure2/",i,"_anno_risk.Rds"))
  dat <- expr_list[[i]][,rownames(anno)]
  print(table(rownames(t(dat)) == rownames(anno)))
  dat_cor <- corAndPvalue(t(dat), anno$RiskScore)
  tmp <- as.data.frame(do.call(cbind, dat_cor))
  colnames(tmp) <- names(dat_cor)
  
  library(GSEABase)
  library(clusterProfiler)
  geneList=tmp$cor
  names(geneList)=rownames(tmp)
  geneList=sort(geneList,decreasing = T)
  d='../../newmodel/Msigdb/'
  gmts <- list.files(d,pattern = 'h')
  gmts
  geneset <- read.gmt(file.path(d,gmts)) 
  gsea_results <- GSEA(geneList, TERM2GENE=geneset, verbose=FALSE)
  
  saveRDS(gsea_results,file = paste0(i,'_GSEA.Rds'))
  
  geneset <- getGmt(file.path(d,gmts)) 
  gsva_results <- GSVA::gsva(as.matrix(dat), geneset, mx.diff=FALSE, verbose=FALSE, method = 'ssgsea')
  
  saveRDS(gsva_results,file = paste0(i,'_GSVA.Rds'))
}

#结果分析与绘图
#GSEA
GSE32537_GSEA <- readRDS("H:/Project/Smart/MMP-1/Figure3/GSE32537_GSEA.Rds")
GSE47460_GSEA <- readRDS("H:/Project/Smart/MMP-1/Figure3/GSE47460_GSEA.Rds")
res32537 <- GSE32537_GSEA@result
res47460 <- GSE47460_GSEA@result
gsea_results <- list(GSE32537 = res32537, GSE47460 = res47460)
pathway_list <- lapply(gsea_results, function(x){
  x <- x$ID
})
lapply(gsea_results, function(x) x[x$ID %in% Reduce(intersect,pathway_list),c('ID','NES')])


dat <- do.call(rbind,gsea_results)
dat$DataSet <- c(rep('GSE32537',nrow(gsea_results$GSE32537)),
                 rep('GSE47460',nrow(gsea_results$GSE47460)))
mycol = ggsci::pal_lancet()(2)
p1 <- ggplot(dat,aes(x=DataSet,y=ID)) #热图绘制
p2 <- p1+scale_fill_gradient2(low = mycol[1], mid = 'white', high = mycol[2],)+
  theme_minimal()+
  geom_point(aes(size=-log10(padjust),
                 fill=NES), shape= 21)+
  theme(axis.text.x =element_text(angle =45,hjust =0.5,vjust = 0.5))+
  xlab(NULL) + ylab(NULL)
p2
ggsave(plot = p2,'GSEA_cor_ipf.pdf', width = 7, height = 6)

#Figure3H
for (i in rownames(res32537)) {
  gseaplot2(GSE32537_GSEA,i,title = i)
  ggsave(paste('GSE32537_',i,'.pdf'),width = 6,height = 5)
}

for (i in rownames(res47460)) {
  gseaplot2(GSE47460_GSEA,i,title = i)
  ggsave(paste('GSE47460_',i,'.pdf'),width = 6,height = 5)
}

#Figure3I
GSE32537_GSVA <- readRDS("H:/Project/Smart/MMP-1/Figure3/GSE32537_GSVA.Rds")
GSE47460_GSVA <- readRDS("H:/Project/Smart/MMP-1/Figure3/GSE47460_GSVA.Rds")
gsva_results <- list(GSE32537_GSVA,GSE47460_GSVA)

library(ggpubr)
for (i in names(gsva_results)) {
  anno <- readRDS(paste0("H:/Project/Smart/MMP-1/Figure2/",i,"_anno_risk.Rds"))
  dat <- gsva_results[[i]][,rownames(anno)]
  for (j in rownames(dat)) {
    print(table(rownames(t(dat)) == rownames(anno)))
    tmp <- data.frame(pathway = dat[j,], RiskScore = anno$RiskScore)
    ggscatter(tmp,y = 'pathway', x = "RiskScore",
              color = "black", size = 3, # 点的颜色，大小
              add = "reg.line",  # 添加回归线
              add.params = list(color = "blue", fill = "lightgray"), # 回归线的调整
              conf.int = TRUE, # 回归线的置信区间
              cor.coef = TRUE, # 添加相关系数
              cor.coeff.args = list(method = "pearson", label.sep = "\n"),
              ylab = j,
              xlab = 'RiskScore')
    ggsave(paste0(i,'_',j,'_cor.pdf'), width = 5, height = 5)
  }
}
