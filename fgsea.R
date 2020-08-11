#BiocManager::install("fgsea")
library(fgsea)
require(data.table)
pathways <- gmtPathways("c5.all.v7.1.symbols.gmt")#msigdb.v7.1.symbols.gmt, ./c2.cp.v7.1.symbols.gmt c5.all.v7.1.symbols.gmt c3.tft.gtrd.v7.1.symbols


stats <- read.delim("mousevsall_statg.txt",row.names=1, dec=',')
sts=stats[,1]
names(sts)=rownames(stats)




fgseaRes <- fgsea(pathways = pathways,stats = sts,minSize  = 2,
                  			maxSize  = 100,nperm=500000)
setorder(fgseaRes,pval)
#setorderv(fgseaRes,'size',-1)

fgseaRes

write.table(as.matrix(fgseaRes), file="./fgseaMousevsAll.txt",
            sep="\t", col.names=T, row.names=T, append = F, quote=FALSE)

