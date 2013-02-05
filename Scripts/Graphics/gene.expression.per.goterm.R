## Load data
args <- commandArgs(TRUE)
#args <- c("LM1_HeLa_all_timep_final.txt", "exp/", "LM1", c("GOsummary_BP_120.csv"))
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
geneExpressionPlot <- function ( gotable, count_table, goname){
	genes=gotable$genes_RefSeq
	goterms=as.list(gotable[1])
	if ( length(genes) > 30 ) { sizeGenes=30 } else { sizeGenes=length(genes) }
	for ( i in 1:sizeGenes){
		splitData=strsplit(as.character(genes[i]),", ")
		length_splitdata=length(splitData[[1]])
		sum=length_splitdata%%2
		gosplit=strsplit(as.character(goname),"_")
		go_ont=gosplit[[1]][2]
		go_timep=strsplit(as.character(gosplit[[1]][3]),"[.]")[[1]][1]
		if(length_splitdata<10){
		png(filename= sprintf("%s%s/%s_%s_%s_%s.png", args[2], args[3], args[3], go_ont, go_timep, goterms[[1]][i])) 
		if (length_splitdata==2){
			par(mfrow=c(2,1), oma = c(1.5,1.5,0,0), mar = c(2.5,2.5,1.5,0.5))
		}
		else if (length_splitdata!=1){
			par(mfrow=c((length_splitdata+sum)/2,2), oma = c(1.5,1.5,0,0), mar = c(2.5,2.5,1.5,0.5))
		}
		for( j in 1:length_splitdata){
				plotExpression(splitData[[1]][j], count_table)
		}
		mtext("Time-points [min]", side=1, outer=TRUE)
		mtext("Nb counts normalized", side=2, outer=TRUE)
		dev.off()
		}
		if(length_splitdata>10){
		length_splitdata_var=length_splitdata/10
		pos=1
		while (length_splitdata_var > 0){
			if(length_splitdata_var > 1){
				png(filename= sprintf("%s%s/%s_%s_%s_%s_%s.png", args[2], args[3], args[3], go_ont, go_timep, goterms[[1]][i], ceiling(length_splitdata_var)))
				par(mfrow=c(5,2), oma = c(1.5,1.5,0,0), mar = c(2.5,2.5,1.5,0.5))
				for( j in pos:(pos+9)){
					plotExpression(splitData[[1]][j], count_table)
					print(j)
				}
				mtext("Time-points [min]", side=1, outer=TRUE)
				mtext("Nb counts normalized", side=2, outer=TRUE)
				dev.off()
				pos=pos+9
				print(pos)
				length_splitdata_var=length_splitdata_var-1
			}
			else{
				sum=length_splitdata_var%%2
				png(filename= sprintf("%s%s/%s_%s_%s_%s_%s.png", args[2], args[3], args[3], go_ont, go_timep, goterms[[1]][i], ceiling(length_splitdata_var)))
				#print((length_splitdata_var*10+sum)/2)
				par(mfrow=c((length_splitdata_var*10+sum)/2,2), oma = c(1.5,1.5,0,0), mar = c(2.5,2.5,1.5,0.5))
				for( j in pos:length_splitdata){
				plotExpression(splitData[[1]][j], count_table)
				}
				print(pos)
				mtext("Time-points [min]", side=1, outer=TRUE)
				mtext("Nb counts normalized", side=2, outer=TRUE)
				dev.off()
				length_splitdata_var=0
			}
		}
	}
	}
}
geneExpressionPlot(gotable, normalized, "GOsummary_BP_120.csv"

for( i in 4:length(args)){
	gotable = read.csv(file = args[i], header=TRUE, row.names=1)
	geneExpressionPlot(gotable, normalized, args[i])
}
