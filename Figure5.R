library(tidyverse)
library(Seurat)
#Figure5AB
load("H:/Project/Smart/MMP/mmp/SingleCell/sce.Rdata")
DimPlot(sce_sub, split.by = 'Diagnosis')
DimPlot(sce_sub, split.by = 'celltype')

#Figure5C
mmp_gene <- grep('^MMP', gene, value = T)
mmp_gene
DotPlot(sce_sub, feature = mmp_gene, cols = c('grey','red'))

#Figure5D
FeaturePlot(sce_sub, feature = mmp_gene, cols = c('grey','red'))

#Figure5E
deg <- FindMarkers(sce_sub,
                   group.by = 'Diagnosis',
                   ident.1 = 'IPF',
                   ident.2 = 'Control',logfc.threshold = 0)

gene = intersect(gene, rownames(deg))
deg$Difference <- deg$pct.1 - deg$pct.2
deg$group <- 'no'
deg[gene,'group'] <- 'yes'
deg$gene_name <- rownames(deg)

library(ggrepel)
ggplot(deg, aes(x = Difference, y = avg_log2FC))+
  geom_point(size = 0.5, aes(color = group))+
  scale_color_manual(values = c('grey','red'))+
  geom_vline(xintercept = 0.0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+
  theme_classic()
