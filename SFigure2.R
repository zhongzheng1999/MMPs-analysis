###SFigure2A
library(missMDA)
library(pheatmap)
load("H:/Project/Smart/MMP/data/GSE47460/GSE47460_expr.Rdata")
metadata <- metadata[metadata$final.diagnosis == 'ILD', ]
lung_function <- metadata[, colnames(metadata)[c(4:8)]]
mean(is.na(lung_function))
lung_function_imputed <- imputePCA(lung_function)$completeObs
pca <- prcomp(lung_function_imputed, scale. = T)

anno <- data.frame(meta_lung_function = -pca$x[,1], gender = as.factor(metadata$Sex), age = metadata$age)
anno$gender[which(anno$gender == '')] <- NA
rownames(anno) <- rownames(lung_function_imputed)
pheatmap(t(lung_function_imputed[order(pca$x[,1]),]),
         cluster_cols = F, annotation_col = anno, scale = 'row',
         show_colnames = F, breaks = seq(-2, 2, length = 101),
         filename = 'GSE47460_meta.pdf', width = 10,height = 6)
saveRDS(anno,'GSE47460.Rds')


mmp_expr <- readRDS("H:/Project/Smart/MMP/mmp/mmp_expr.Rds")
dat <- mmp_expr$GSE32537
meta <- readRDS("H:/Project/Smart/MMP/mmp/GSE32537meta.Rds")
meta <- meta[meta$final.diagnosis == 'IPF',]
colnames(meta)
colnames(meta) <- c('age',"gender","final.diagnosis","DLCO",
                    "FVC","pack.years","preservative",
                    "quit","repository","rin","smoking.status",
                    "george.total.score","tissue.source")

dat <- dat[rownames(meta)]
lung_function <- meta[, colnames(meta)[c(4,5,12)]]
mean(is.na(lung_function))
lung_function_imputed <- imputePCA(lung_function)$completeObs
pca <- prcomp(lung_function_imputed, scale. = T)

anno <- data.frame(meta_lung_function = pca$x[,1], gender = as.factor(meta$gender), age = meta$age)
anno$gender[which(anno$gender == '')] <- NA
rownames(anno) <- rownames(lung_function_imputed)
pheatmap(t(lung_function_imputed[order(pca$x[,1]),]),
         cluster_cols = F, annotation_col = anno, scale = 'row',
         show_colnames = F, breaks = seq(-2, 2, length = 101),
         filename = 'GSE32537_meta.pdf', width = 9,height = 2.5)
saveRDS(anno,'GSE32537.Rds')

#SFigure 2C
library(ggpubr)
dat <- t(ipf_expr$GSE47460)
rownames(anno) == rownames(dat)
for (i in rownames(ipf_expr$GSE47460)) {
  group <- ifelse(dat[,i] > median(dat[,i]), 'High','Low')
  tmp <- data.frame(meta = anno$meta_lung_function,
                    exp = dat[,i],
                    group = group)
  ggviolin(tmp, x = 'group', y = 'meta', fill = 'group', palette = 'lancet',
           add = 'boxplot', add.params = list(fill = 'white'))+
    stat_compare_means()
  ggsave(paste0(i,'_function_vlinplot.pdf'),width = 4, height = 5)
  
}


#SFigure2D,E
lasso_fea <- readRDS("H:/Project/Smart/MMP-1/Figure2/lasso_fea.Rds")
train <- as.data.frame(dat[,lasso_fea]) 
dim(train)
train[1:4,1:4]
x <- scale(train)
rownames(train) == rownames(anno)
train$group <- anno$meta_lung_function
# (y <- ifelse(train$group > 0, 1, 0)) #把分组信息换成01
(y <- train$group)
fit = glmnet(x, y, family = "gaussian", alpha = 0, lambda = NULL)
pdf("A_rid.pdf", width = 5, height = 5)
plot(fit, xvar = "dev", label = TRUE)
cvfit = cv.glmnet(x, y, alpha = 0,
                  nfold=10, #例文描述：10-fold cross-validation
                  family = "gaussian", type.measure = "class")
plot(cvfit)

dev.off()

cvfit$lambda.min
myCoefs <- coef(cvfit, s="lambda.min");myCoefs


#SFigure2F,G

anno <- readRDS("H:/Project/Smart/MMP-1/Figure2/GSE47460_anno_risk.Rds")
rownames(metadata) == rownames(anno)
ggboxplot(anno, x = 'gender', y = 'RiskScore', fill = 'gender', palette = 'lancet')+
  stat_compare_means()
ggsave('GSE47460_Sex_risk.pdf',width = 5, height = 5)

ggscatter(anno,x = 'age', y = "RiskScore",
          color = "black", size = 3, # 点的颜色，大小
          add = "reg.line",  # 添加回归线
          add.params = list(color = "blue", fill = "lightgray"), # 回归线的调整
          conf.int = TRUE, # 回归线的置信区间
          cor.coef = TRUE, # 添加相关系数
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          xlab = 'Meta Lung Function',
          ylab = 'Predictive Lung Function')
ggsave('GSE47460_age_predictive_cor.pdf', width = 5, height = 5)