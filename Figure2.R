#Figure2B
ipf_expr <- readRDS("H:/Project/Smart/MMP-1/ipf_expr.Rds")
colnames(ipf_expr$GSE47460) == rownames(anno)
dat <- ipf_expr$GSE47460
res <- t(apply(dat, 1, function(x){
  afit <- lm(x ~ anno$gender + anno$age + anno$meta_lung_function) 
  tmp <- coefficients(summary(afit))
  tmp[nrow(tmp),]
}))
colnames(res) <- c('estimate', 'std_error', 't_value', 'p_value')
res <- as.data.frame(res)
res$Gene.name <- rownames(res)
res$adjusted_pval <- p.adjust(res$p_value, method = 'BH')
write.csv(res,'GSE47460_gene_lungfuntion.csv')

sig <- which(res$adjusted_pval < 0.25)

library(ggrepel)
aframe <- data.frame(coef = res$estimate, pval = -log10(res$p_value), gene = res$Gene.name)
ggplot(aframe, aes(x = coef, y = pval, color = pval)) + geom_point(alpha = 0.5) +
  xlab('Coefficient') + ylab('Significance [-log10 p-value]') +
  scale_colour_gradient2(low = "black", mid = 'grey', high = "red") + 
  geom_text_repel(data = aframe[sig,], aes(label = gene))
ggsave('coef_pvalue.pdf',width = 4, height = 5)


#Figure2C,D
library(tidyverse)
library(glmnet)
library(VennDiagram)
library(e1071)
library(caret)
library(randomForest)
ipf_expr <- readRDS("H:/Project/Smart/MMP-1/ipf_expr.Rds")
dat <- t(ipf_expr$GSE47460)
all_mmpsig <- c('MMP9','MMP7','MMP25','MMP28','MMP19','MMP17','MMP14','MMP11','MMP10','MMP1')
train <- as.data.frame(dat[,all_mmpsig])

dim(train)

train[1:4,1:4]
x <- scale(train)
rownames(train) == rownames(anno)
train$group <- anno$meta_lung_function

# (y <- ifelse(train$group > 0, 1, 0)) #把分组信息换成01
(y <- train$group)
#library(glmnet)
fit = glmnet(x, y, family = "gaussian", alpha = 1, lambda = NULL)

# 画A图
pdf("A_lasso.pdf", width = 5, height = 5)
plot(fit, xvar = "dev", label = TRUE)

cvfit = cv.glmnet(x, y, alpha = 1,
                  nfold=10, #例文描述：10-fold cross-validation
                  family = "gaussian", type.measure = "class")
plot(cvfit)

dev.off()

cvfit$lambda.min

# 获取LASSO选出来的特征
myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
(lasso_fea <- lasso_fea[-1])

write.csv(lasso_fea,"feature_lasso.csv")
saveRDS(lasso_fea,"lasso_fea.Rds")
save(fit,cvfit,file = 'lasso_fit.Rdata')

#Figure2E
fit = glmnet(x, y, family = "gaussian", alpha = 0, lambda = cvfit$lambda.min)
coef_dat <- as.data.frame(as.matrix(coef(fit)))
coef_dat$Variable <- rownames(coef_dat)
coef_dat <- coef_dat[-1,]
coef_dat$type <- 'up'
coef_dat$type[3] <- 'down'
ggbarplot(coef_dat, y = 's0', x = 'Variable', fill = 'type', palette = 'lancet',
          sort.val = 'asc', rotate = T, ylab = 'Coeffecients')
ggsave('Coeddecients.pdf',width = 6, height = 7)

#Figure2F

predCV <- predict(fit, newx = x,
                  # s = "lambda.min",
                  type = "response")
rownames(predCV) == rownames(x)
df <- data.frame(group = y, x, Combined = predCV[,1])

ggscatter(df,x = 'group', y = "Combined",
          color = "black", size = 3, # 点的颜色，大小
          add = "reg.line",  # 添加回归线
          add.params = list(color = "blue", fill = "lightgray"), # 回归线的调整
          conf.int = TRUE, # 回归线的置信区间
          cor.coef = TRUE, # 添加相关系数
          cor.coeff.args = list(method = "pearson", label.sep = "\n"),
          xlab = 'Meta Lung Function',
          ylab = 'Predictive Lung Function')
ggsave('meta_predictive.pdf', width = 5, height = 5)
