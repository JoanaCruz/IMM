# LM_rep1=read.table("DESeq_rep1_Control_ALL.csv", sep=",", header=T)
# LM_rep2=read.table("DESeq_rep2_Control_ALL.csv", sep=",", header=T)
# LM_rep3=read.table("DESeq_rep3_Control_ALL.csv", sep=",", header=T)
# LM_rep4=read.table("DESeq_rep4_Control_ALL.csv", sep=",", header=T)


#LM_rep1_rep2=read.table("DESeq_rep1_rep2_ALL.csv", sep=",", header=T)
args <- commandArgs(TRUE)

plotDE <- function( res ){
    plot( 
	res$baseMean, 
	res$log2FoldChange,
	log="x", pch=20, cex=.3, 
	col = ifelse( res$padj < .1, "red", "black" ),
	ylim = c(-10,10)
#	xlim = c(0,exp(7))
	)
}

# png("dispersion.plot.rep3.png", width=940,height=840)
# 
# par(mfrow=c(4,1), mar=c(5,5,2,2))
# plotDE(LM_rep1)                                                                     
# abline(h=0, col="red3", lwd=2)
# title("rep3_1 vs control")
# plotDE(LM_rep2)
# abline(h=0, col="red3", lwd=2)
# title("rep3_2 vs control")
# plotDE(LM_rep3)
# abline(h=0, col="red3", lwd=2)
# title("rep3_3 vs control")
# plotDE(LM_rep4)
# abline(h=0, col="red3", lwd=2)
# title("rep3_4 vs control")
# dev.off()


files = list.files(args[1], pattern="*_ALL.csv")

for ( i in 1:length(files) ) {
	name=gsub("DESeq_","",files[i])
	name1=gsub("_ALL.csv","",name)
	data=read.table(sprintf("%s%s",args[1],files[i]), sep=",", header=T)
	png(sprintf("%sDispersion_plot_%s.png", args[2], name1), width=540, height=440)
	plotDE(data)                                                                     
	abline(h=0, col="red3", lwd=2)
	title(sprintf("%s",name1))
	dev.off()
}