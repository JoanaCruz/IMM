args <- commandArgs(TRUE)

#load table with counts
count_table = read.table(args[1], header=T, sep="\t", row.names=1)

analysis=c("rep1.vs.rep2", "public_rep.vs.control", "rep1.vs.control", "rep2.vs.control")
library( "DESeq" )

for ( cond in analysis){

	# define data design & table with counts
	if (cond=="rep1.vs.rep2"){
		dataDesign = data.frame(
			row.names = c( colnames( count_table )[2], colnames( count_table )[3]),
			condition = c( "rep1", "rep2"),
			libType = rep("paired-end", 2))
		count_table_analysis = count_table[,c( colnames( count_table )[2], colnames( count_table )[3])]
		methodology='blind'
		sharingM='fit-only'
		fitT='local'
	}
	else if (cond=="public_rep.vs.control"){
		dataDesign = data.frame(
			row.names = colnames( count_table ),
			condition = c( "Control_HeLa", rep("rep",2) ),
			libType = rep("paired-end",3))
		count_table_analysis = count_table[,c( colnames( count_table )[1], colnames( count_table )[2], colnames( count_table )[3])]
		methodology='pooled'
		sharingM='maximum'
		fitT='parametric'
	}
	else if (cond=="rep1.vs.control"){
		dataDesign = data.frame(
			row.names = c( colnames( count_table )[1], colnames( count_table )[2]),
			condition = c( "Control", "rep1"),
			libType = rep("paired-end",2))
		count_table_analysis = count_table[,c( colnames( count_table )[1], colnames( count_table )[2])]
		methodology='blind'
		sharingM='fit-only'
		fitT='local'
	}
	else if (cond=="rep2.vs.control"){
		dataDesign = data.frame(
			row.names = c( colnames( count_table )[1], colnames( count_table )[3]),
			condition = c( "Control", "rep2"),
			libType = rep("paired-end",2))
		count_table_analysis = count_table[,c( colnames( count_table )[1], colnames( count_table )[3])]
		methodology='blind'
		sharingM='fit-only'
		fitT='local'
	}
	
	conditions=dataDesign$condition
	libraryType=dataDesign$libType
		
	
	data = newCountDataSet(count_table_analysis, conditions)


# 1. NORMALIZATION
# estimate the size factors (which is a representation of the differences in coverage between replicates)
# The effective library size information is called the size factors vector, since the package only needs to know the relative library sizesdata
	data = estimateSizeFactors(data)

#write(sizeFactors(data),file=sprintf("LM_HeLa_%s.txt",args[2]))

# 2. VARIANCE

# The inference in DESeq relies on an estimation of the typical relationship between the data’s variance and their mean, or, equivalently, between the data’s dispersion and their mean.

# Verification of the dispersion is an aspect of data quality control: how well do the data accord to the expectations of our analytic approach?

# estimateDispersions: First, it estimates a dispersion value for each gene, then it fits a curve through the estimates. Finally, it assigns to each gene a dispersion value, using a choice between the per-gene estimate and the fitted value.

	data = estimateDispersions( data, method=methodology, sharingMode=sharingM, fitType=fitT)

# Now, we can attempt to find differential expression:
	res = nbinomTest( data, unique(conditions)[1], unique(conditions)[2] )

	resp=subset(res,res$padj<0.1)

	res_pvalue=resp[ order(resp$pval), ]

	res_down=resp[ order( resp$foldChange, -resp$baseMean ), ]

	res_up=resp[ order( -resp$foldChange, -resp$baseMean ), ]
	
	write.csv(res_pvalue, file = sprintf("%sDESeq_%s_pvalue.csv", args[2], cond))

	write.csv(res_up, file = sprintf("%sDESeq_%s_up.csv", args[2], cond))

	write.csv(res_down, file = sprintf("%sDESeq_%s_down.csv", args[2], cond))

	write.csv(res, file = sprintf("%sDESeq_%s_ALL.csv", args[2], cond))
}