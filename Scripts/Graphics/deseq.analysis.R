LM_20 = read.table("DESeq_LM_HeLa_20_ALL.txt")
LM_60 = read.table("DESeq_LM_HeLa_60_ALL.txt")
LM_120 = read.table("DESeq_LM_HeLa_120_ALL.txt")
LM_240 = read.table("DESeq_LM_HeLa_240_ALL.txt")


plotDE <- function( res )
   plot( 
      res$baseMean, 
      res$log2FoldChange,
      log="x", pch=20, cex=.3, 
      col = ifelse( res$padj < .1, "red", "black" ) )


par(mfrow=c(2,2))
plotDE(LM_20)
title(main="LM_20")
plotDE(LM_60)
title(main="LM_60")
plotDE(LM_120)
title(main="LM_120")
plotDE(LM_240)
title(main="LM_240")
#par(mfrow=c(4,1))

#hist(LM_20$padj, breaks=100, col="skyblue", border="slateblue", main="")
#hist(LM_60$padj, breaks=100, col="skyblue", border="slateblue", main="")
#hist(LM_120$padj, breaks=100, col="skyblue", border="slateblue", main="")
#hist(LM_240$padj, breaks=100, col="skyblue", border="slateblue", main="")


#addmargins( table( LM_20 = LM_20$padj < .1, LM_60 = LM_60$padj < .1 ) )


#pval <- rbind(table(LM_20$padj < .1), table(LM_60$padj < .1), table(LM_120$padj < .1), table(LM_240$padj< .1))

#rownames(pval) <- c("LM_20","LM_60","LM_120","LM_240")

#write.table(pval, file = "p_va_refseq.txt")


#LM20_p=subset(LM_20,LM_20$padj<0.1)

#LM20_pvalue=LM20_p[ order(LM20_p$pval), ]

#write.table(LM20_pvalue, file = "DESeq_LM_HeLa_20_pvalue.txt")

