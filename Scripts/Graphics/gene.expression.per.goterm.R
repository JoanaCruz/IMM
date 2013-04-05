## Load data
args <- commandArgs(TRUE)

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
	plot(y,time_counts, xlab="Time [minutes]", ylab="Normalized expression level (Nb of counts)", main=gene, xlim=c(0,245), pch=19, col="blue")
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
				dev.off()
			}
		}
	}
}

if(length(args)>3){
	# Plots only the differentially expressed genes associated with the go terms found by GOStats
	for( i in 3:length(args)){
		gotable = read.csv(file = args[i], header=TRUE, row.names=1)
		geneExpressionPlot(gotable, normalized)
	}
} else {
 	#Plots the gene expression profile associated with all genes in the analysis
 	mrna_data=file(args[2], 'r')
 	line=readLines(mrna_data, warn=FALSE)
 	genes=strsplit(line, " ")
	for ( i in 1:length(genes[[1]]) ){
		gene=genes[[1]][i]
		sample = strsplit(strsplit(as.character(args[1]),"/")[[1]][3],"_")[[1]][1]
		png(filename= sprintf("%s%s_%s.png", args[3], sample, gene))
		plotExpression(gene, normalized)
		dev.off()
	}
}
# 	for( pos in 1:nrow(normalized) ){
# 		y = c(0,20,60,120,240)
# 		counts = as.matrix(normalized[pos,])
# 		gene = row.names(normalized)[pos]
# 		
# 		
# 		# Plot data
# 		
# 		plot(y,counts, xlab="Time [minutes]", ylab="Normalized expression level (Nb of counts)", main=gene, xlim=c(0,245), pch=19, col="blue")
# 		dev.off()
# 	}
#}
	
