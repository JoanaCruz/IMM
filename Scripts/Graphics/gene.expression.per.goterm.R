## Load data
args <- commandArgs(TRUE)
count_tabe = read.table(file = args[1], header=TRUE, row.names=1)
gotable = read.csv(file = args[2], header=TRUE, row.names=1)

library( "DESeq" )


## Define data design
dataDesign = data.frame(
        row.names = colnames( count_table ),
        condition = colnames( count_table ),
        libType = c("paired-end","paired-end","paired-end","paired-end","paired-end"))

conditions=dataDesign$condition

libraryType=dataDesign$libType

data = newCountDataSet(count_table, conditions)


## 1. NORMALIZATION
# estimate size factors (which is a representation of the differences in coverage between replicates)
# The effective library size information is called the size factors vector, since the package only needs to know the relative library sizesdata

data = estimateSizeFactors(data)

normalized = counts(data, normalized=TRUE)


## Gene expression graphic

#Gene profile for a single gene in table of counts
plotExpression <- function ( gene, table ){
	pos = row.names(table) == gene
	time_counts = as.matrix(table[pos,])
	y = c(0,20,60,120,240)
	plot(y,time_counts, xlab=NA, ylab=NA, main=gene, xlim=c(0,250), pch=19, col="blue")
	lines(y,time_counts)
}


#Plots the expression profile of the genes associated with each one of the go terms in gotable
geneExpressionPlot <- function ( gotable, count_table){
	genes=gotable$genes_RefSeq
	goterms=gotable$GOBPID
	for ( i in i:length(genes)){ 
		splitData=strsplit(as.character(genes[i]),", ")
		length_splitdata=length(splitData[[1]])
		sum=length_splitdata%%2
		if (length_splitdata!=1){
			par(mfrow=c((length_splitdata+sum)/2,2), oma = c(1.5,1.5,0,0), mar = c(2.5,2.5,1.5,0.5))
		}
		png(filename= sprintf("%s%s.png", arg[3], goterms[i])) 
		for( j in 1:length_splitdata){
				plotExpression(splitData[[1]][j], count_table)
		}
		mtext("Time-points [min]", side=1, outer=TRUE)
		mtext("Nb counts normalized", side=2, outer=TRUE)
		dev.off()
	}
}

