#!/usr/bin/env Rscript
library("edgeR")
require(data.table)

# данные : data = data.table, колонки: region, h_liver, h_brain, m_liver, m_brain
data=fread("comparison.txt")
#data <- na.omit(data)

y <- DGEList(counts= data.matrix(data[,c(2:33)]), genes=data[,1])
rownames(y$counts) <- rownames(y$genes) <- y$genes$gene

#filter
keep <- rowSums(cpm(y)>2) >= 4 #set cpm
y <- y[keep, ,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)

organism <- factor(rep(c('h','m'), c(15,17)))
tissue <- factor(substr(colnames(y),1,nchar(colnames(y))-2))
organism <- relevel(organism, ref="h")
tissue <- relevel(tissue, ref="astrocite_cerebell")
group = tissue:organism

design <- model.matrix(~0+group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
cntrs=c(0, 0,-1,0,0,0,1,0,0,0)

qlf <- glmQLFTest(fit,contrast=cntrs)
all_genes_res <- as.data.table(topTags(qlf, n = nrow(y), sort.by = "none"))

write.table(all_genes_res, file="~/diffexp.txt",
            sep="\t", col.names=T, row.names=T, append = F, quote=FALSE)