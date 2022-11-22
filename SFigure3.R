library(survival)
library(survminer)
library(patchwork)
load("H:/Project/Smart/MMP/data/BAL_IPF_cohort.Rdata")
lasso_fea <- readRDS("H:/Project/Smart/MMP-1/Figure2/lasso_fea.Rds")
#全局变量
methon <- 'best' #median mean best
GSE_num <- 'all' #"LEUVEN"   "SIENA"    "Freiburg"


for (gene in lasso_fea) {
  #构建数据框
  dat <- data.frame(t(exp[[GSE_num]][gene,]),
                    meta[[GSE_num]]['State'],
                    meta[[GSE_num]]['Days'])
  dat$State <- as.numeric(dat$State)
  dat$Days <- as.numeric(dat$Days)
  
  
  #正式流程
  #判断methon
  res.cut <- surv_cutpoint(dat,
                           time ="Days",
                           event = "State",
                           variables = gene,
                           minprop = 0.3)
  risk <- surv_categorize(res.cut)
  
  #生存曲线
  dat$risk <- risk[,gene]
  surv_object = Surv(dat$Days, dat$State)
  fit1 <- survfit(surv_object ~ risk, data = dat)
  summary(fit1)
  
  p <- ggsurvplot(fit1,palette = c("#E7B800", "#2E9FDF"),
                  risk.table =TRUE,pval =TRUE,xlab ="Time in Days",
                  ggtheme =theme_light(),title = gene,
                  ncensor.plot = F)
  p1 <- p$plot/p$table + plot_layout(heights = c(3, 1))
  ggsave(plot = p1, filename = paste0(GSE_num, '_', gene, '_Survival.pdf'), width = 6, height = 5)
  
}


#单因素COX
load("H:/Project/Smart/MMP/data/BAL_IPF_cohort.Rdata")

meta$all <- rbind(meta$LEUVEN,meta$Freiburg,meta$SIENA)
meta$all$group <- c(rep('LEUVEN',nrow(meta$LEUVEN)),
                    rep('Freiburg',nrow(meta$Freiburg)),
                    rep('SIENA',nrow(meta$SIENA)))
exp$all_un <- cbind(exp$LEUVEN,exp$Freiburg,exp$SIENA)
colnames(exp$all_un) == rownames(meta$all)
library(sva)
exp$all <- ComBat(dat = exp$all_un, batch = meta$all$group)
exp$all <- as.data.frame(exp$all)

GSE_num <- 'all' #"LEUVEN"   "SIENA"    "Freiburg"
#构建数据框
dat <- data.frame(t(exp[[GSE_num]][lasso_fea,]),
                  meta[[GSE_num]]['State'],
                  meta[[GSE_num]]['Days'])
dat$State <- as.numeric(dat$State)
dat$Days <- as.numeric(dat$Days)
res.cox <- coxph(Surv(Days, State) ~ MMP19, data = dat)
res.cox
#对所有特征做单因素cox回归分析

#分别对每一个变量，构建生存分析的公式
covariates <- colnames(dat)[-c(9,10)]

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Days, State)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = dat)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$coefficients[,5], digits=2)
                         #获取HR
                         HR <-signif(x$conf.int[,1], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,3], 2)
                         HR.confint.upper <- signif(x$conf.int[,4], 2)
                         HR_all <- paste0(HR, " (", 
                                          HR.confint.lower, "-", HR.confint.upper, ")")
                         cof <- signif(x$coef[1], 2)
                         res<-c(cof,p.value,HR_all,HR,HR.confint.lower,HR.confint.upper)
                         names(res)<-c("coef","p.value","HR (95% CI for HR)","HR","HR.confint.lower","HR.confint.upper")
                         return(res)
                       })
#转换成数据框，并转置
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
write.table(file=paste0(GSE_num,"_univariate_cox_result.txt"),res,quote=F,sep="\t")

#作图森林图
library(grid)
library(magrittr)
library(checkmate)
library(forestplot)

res$coef <- round(as.numeric(res$coef),3)
res$p.value <- as.numeric(res$p.value)
res$HR <- as.numeric(res$HR)
res$HR.confint.lower <- as.numeric(res$HR.confint.lower)
res$HR.confint.upper <- as.numeric(res$HR.confint.upper)

mydata <- res
mydata$Gene <- rownames(mydata)
# mydata <- mydata[mydata$p.value < 0.05,]
# write.csv(mydata,'univariate_cox_p0.05.csv')
mydata <- mydata[order(mydata$p.value),]
mydata <- mydata[c(7,1:6)]


pdf(paste0(GSE_num,'_forestplot.pdf'),width = 5,height = 3.5, onefile=FALSE)
forestplot(labeltext=as.matrix(mydata[,1:4]),#只展示前面的三列
           mean= mydata$HR,#HR值
           lower= mydata$HR.confint.lower,#95%CI下限
           upper= mydata$HR.confint.upper,#95%CI上限
           zero=1,#HR值的位置，如果是线性回归则写0
           boxsize=0.1,#中间方框的大小
           xticks=seq(-1,5,1),#x轴刻度
           lwd.zero=2,#中间竖线的宽度
           lwd.ci=2, 
           col=fpColors(box='orange',lines = 'orange',zero = 'gray'),#颜色，box，lines和zero分别是方框，线条，中间竖线的颜色
           xlab="HR",#x轴标签
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),xlab  = gpar(cex = 0.8),
                            cex = 0.9),#设置字体大小
           lty.ci = "solid",
           title = GSE_num, #标题
           graph.pos=2#中间竖线的位置
)
dev.off()


#模型构建
library(glmnet)
library(survival)
GSE_num <- 'all' #"LEUVEN"   "SIENA"    "Freiburg"
set.seed(100)

unicox <- read.table('all_univariate_cox_result.txt',header = T,sep = '\t')
unicox <- unicox[unicox$p.value < 0.05,]
gene <- rownames(unicox)

dat <- data.frame(t(exp[[GSE_num]][gene,]),
                  meta[[GSE_num]]['State'],
                  meta[[GSE_num]]['Days'])
dat$Days <- as.numeric(dat$Days)
dat$State <- as.numeric(dat$State)

y <- data.matrix(Surv(time = dat$Days,event = dat$State))


marker.exp <- dat[,gene]
#构建模型
fit <- glmnet(x = as.matrix(marker.exp),y,family = 'cox',alpha = 0)
pdf(paste0(GSE_num, '_ridge1.pdf'))
plot(fit,xvar = 'lambda',label = T)
dev.off()

ridge_fit <- cv.glmnet(x = as.matrix(marker.exp),y,family = 'cox',type.measure = 'deviance',alpha = 0)
pdf(paste0(GSE_num, '_ridge2.pdf'))
plot(ridge_fit,label = T)
dev.off()
save(fit,ridge_fit,file = paste0(GSE_num,'_ridge_gene.Rdata'))


fit <- glmnet(x = as.matrix(marker.exp),y,family = 'cox',alpha = 0, lambda = ridge_fit$lambda.min)
#saveRDS(fit,'fit_model.Rds')
dat$score <- as.numeric(predict(fit,newx = as.matrix(marker.exp)))


res.cut <- surv_cutpoint(dat, 
                         time ="Days", 
                         event = "State", 
                         variables = "score", 
                         minprop = 0.3)  
risk <- surv_categorize(res.cut)
dat$risk <- risk$score

surv_object = Surv(dat$Days, dat$State)
fit1 <- survfit(surv_object ~ risk, data = dat)
summary(fit1)
pdf('merge_survive.pdf',onefile = F, width = 6, height = 5)
ggsurvplot(fit1,palette = c("#E7B800", "#2E9FDF"),
           risk.table =TRUE,pval =TRUE,
           conf.int =TRUE,xlab ="Time in months", 
           ggtheme =theme_light(), 
           ncensor.plot = F)
dev.off()

#单因素
rownames(dat) == rownames(meta$all)
dat$Age <- as.numeric(meta$all$`age:ch1`)
dat$Sex <- as.numeric(meta$all$`sex (0=female, 1=male):ch1`)
covariates <- colnames(dat)[-c(1:7,9)]

univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(Days, State)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = dat)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$coefficients[,5], digits=2)
                         #获取HR
                         HR <-signif(x$conf.int[,1], digits=2);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,3], 2)
                         HR.confint.upper <- signif(x$conf.int[,4], 2)
                         HR_all <- paste0(HR, " (", 
                                          HR.confint.lower, "-", HR.confint.upper, ")")
                         cof <- signif(x$coef[1], 2)
                         res<-c(cof,p.value,HR_all,HR,HR.confint.lower,HR.confint.upper)
                         names(res)<-c("coef","p.value","HR (95% CI for HR)","HR","HR.confint.lower","HR.confint.upper")
                         return(res)
                       })
#转换成数据框，并转置
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
write.table(file=paste0(GSE_num,"score_univariate_cox_result.txt"),res,quote=F,sep="\t")

#作图森林图
library(grid)
library(magrittr)
library(checkmate)
library(forestplot)

res$coef <- round(as.numeric(res$coef),3)
res$p.value <- as.numeric(res$p.value)
res$HR <- as.numeric(res$HR)
res$HR.confint.lower <- as.numeric(res$HR.confint.lower)
res$HR.confint.upper <- as.numeric(res$HR.confint.upper)

mydata <- res
mydata$Gene <- rownames(mydata)
# mydata <- mydata[mydata$p.value < 0.05,]
# write.csv(mydata,'univariate_cox_p0.05.csv')
mydata <- mydata[order(mydata$p.value),]
mydata <- mydata[c(3,2,1),c(7,1:6)]
rownames(mydata)[3] <- 'Score'

pdf(paste0(GSE_num,'score_forestplot.pdf'),width = 5,height = 3.5, onefile=FALSE)
forestplot(labeltext=as.matrix(mydata[,1:4]),#只展示前面的三列
           mean= mydata$HR,#HR值
           lower= mydata$HR.confint.lower,#95%CI下限
           upper= mydata$HR.confint.upper,#95%CI上限
           zero=1,#HR值的位置，如果是线性回归则写0
           boxsize=0.1,#中间方框的大小
           xticks=seq(-1,5,1),#x轴刻度
           lwd.zero=2,#中间竖线的宽度
           lwd.ci=2, 
           col=fpColors(box='orange',lines = 'orange',zero = 'gray'),#颜色，box，lines和zero分别是方框，线条，中间竖线的颜色
           xlab="HR",#x轴标签
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),xlab  = gpar(cex = 0.8),
                            cex = 0.9),#设置字体大小
           lty.ci = "solid",
           title = GSE_num, #标题
           graph.pos=2#中间竖线的位置
)
dev.off()

fit1 <- coxph(Surv(Days, State)~ Age+Sex+score, data = dat)

x <- summary(fit1)
#获取p值
p.value<-signif(x$coefficients[,5], digits=2)
#获取HR
HR <-signif(x$conf.int[,1], digits=2);
#获取95%置信区间
HR.confint.lower <- signif(x$conf.int[,3], 2)
HR.confint.upper <- signif(x$conf.int[,4], 2)
HR_all <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
cof <- signif(x$coef[,1], 2)
res<-data.frame(cof,p.value,HR_all,HR,HR.confint.lower,HR.confint.upper)
names(res)<-c("coef","p.value","HR (95% CI for HR)","HR","HR.confint.lower","HR.confint.upper")

res$coef <- round(as.numeric(res$coef),3)
res$p.value <- as.numeric(res$p.value)
res$HR <- as.numeric(res$HR)
res$HR.confint.lower <- as.numeric(res$HR.confint.lower)
res$HR.confint.upper <- as.numeric(res$HR.confint.upper)

mydata <- res
mydata$Gene <- rownames(mydata)
# mydata <- mydata[mydata$p.value < 0.05,]
# write.csv(mydata,'univariate_cox_p0.05.csv')
mydata <- mydata[,c(7,1:6)]
rownames(mydata)[3] <- 'Score'

pdf(paste0(GSE_num,'_muiltscore_forestplot.pdf'),width = 5,height = 3.5, onefile=FALSE)
forestplot(labeltext=as.matrix(mydata[,1:4]),#只展示前面的三列
           mean= mydata$HR,#HR值
           lower= mydata$HR.confint.lower,#95%CI下限
           upper= mydata$HR.confint.upper,#95%CI上限
           zero=1,#HR值的位置，如果是线性回归则写0
           boxsize=0.1,#中间方框的大小
           xticks=seq(-1,5,1),#x轴刻度
           lwd.zero=2,#中间竖线的宽度
           lwd.ci=2, 
           col=fpColors(box='orange',lines = 'orange',zero = 'gray'),#颜色，box，lines和zero分别是方框，线条，中间竖线的颜色
           xlab="HR",#x轴标签
           lwd.xaxis =1,
           txt_gp = fpTxtGp(ticks = gpar(cex = 0.85),xlab  = gpar(cex = 0.8),
                            cex = 0.9),#设置字体大小
           lty.ci = "solid",
           title = GSE_num, #标题
           graph.pos=2#中间竖线的位置
)
dev.off()


library(timeROC)
result <-with(dat, timeROC(T=Days,
                           delta=State,
                           marker=score,
                           cause=1,
                           times=c(365,730,1095),
                           iid = TRUE))
#identical(c(result$TP[,1],result$TP[,2],result$TP[,3]),as.numeric(result$TP))
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(365,730,1095)),each = nrow(result$TP)))

library(ggplot2)

pdf(paste0(GSE_num, '_timeroc.pdf'),onefile = F)
ggplot() + 
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5"),
                     labels = paste0("AUC of ",c(1,2,3),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
dev.off()