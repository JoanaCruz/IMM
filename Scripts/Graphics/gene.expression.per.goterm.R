## Load data
args <- commandArgs(TRUE)
#args <- c("LM1_HeLa_all_timep.txt", "exp/", "LM1", c("GOsummary_BP_120.csv"))
count_table = read.table(file = args[1], header=TRUE, row.names=1)


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
	plot(y,time_counts, xlab=NA, ylab=NA, main=gene, xlim=c(0,245), pch=19, col="blue")
	lines(y,time_counts)
}


#Plots the expression profile of the genes associated with each one of the go terms in gotable
geneExpressionPlot <- function ( gotable, count_table){
	all_genes=gotable$genes_RefSeq
	goterms=as.list(gotable[1])
	sizeGenes=length(all_genes)
	for ( i in 1:sizeGenes){
		splitData=strsplit(as.character(all_genes[i]),", ")
		length_splitdata=length(splitData[[1]])
		for( j in 1:length_splitdata){
			gene=splitData[[1]][j]
			files=list.files(args[2])
			if( !is.element(sprintf("%s.png",gene), files) ){
				png(filename= sprintf("%s%s.png", args[2], gene)) 
				plotExpression(gene, count_table)
				mtext("Time-points [min]", side=1, outer=TRUE)
				mtext("Nb counts normalized", side=2, outer=TRUE)
				dev.off()
			}
		}
	}
}

for( i in 3:length(args)){
	gotable = read.csv(file = args[i], header=TRUE, row.names=1)
	geneExpressionPlot(gotable, normalized)
}
