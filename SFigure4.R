gene <- openxlsx::read.xlsx('1-s2.0-S1074761318301213-mmc7.xlsx',sheet = 1)
gene <- gene[1:78,]

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

dat_imm <- dat[gene$HGNC.Symbol,Sample]
dat_imm <- na.omit(dat_imm)
metadata <- metadata[Sample,]
rownames(metadata) == rownames(GSE47460_anno_risk)
metadata$RiskScore = GSE47460_anno_risk$RiskScore

n <- t(scale(t(dat_imm)))
n[n < -2] = -2
n[n > 2] = 2
anno_row <- data.frame(row.names = gene$HGNC.Symbol, gene[5:7])
anno_row <- anno_row[rownames(n),]
pheatmap::pheatmap(n[,order(metadata$RiskScore)],
                   annotation_col = metadata,annotation_row = anno_row,
                   cluster_cols = T, show_colnames = F,
                   filename = 'immgene_pheatmap_row_ipf.pdf',
                   width = 10, height = 11)


gene <- readRDS("H:/Project/Smart/MMP-1/Figure2/lasso_fea.Rds")
dat <- as.data.frame(t(dat[Sample]))[gene]
rownames(dat) == rownames(GSE47460_anno_risk)
dat$RiskScore <- GSE47460_anno_risk$RiskScore
xcell <- t(dat_imm)
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
ggsave(plot = p2,'immgene_cor_ipf.pdf', width = 6, height = 12)
