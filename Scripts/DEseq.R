args <- commandArgs(TRUE)

#load table with counts
count_table = read.table(args[1], header=T, sep="\t", row.names=1)

# define data design
dataDesign = data.frame(
        row.names = colnames( count_table ),
        condition = c( colnames( count_table )[1], colnames( count_table )[2]),
        libType = c("paired-end","paired-end"))

conditions=dataDesign$condition

libraryType=dataDesign$libType

library( "DESeq" )

data = newCountDataSet(count_table, conditions)


# 1. NORMALIZATION
# estimate the size factors (which is a representation of the differences in coverage between replicates)
# The effective library size information is called the size factors vector, since the package only needs to know the relative library sizesdata

data = estimateSizeFactors(data)

write(sizeFactors(data),file=sprintf("LM_HeLa_%s.txt",args[2]))

# 2. VARIANCE

# The inference in DESeq relies on an estimation of the typical relationship between the data’s variance and their mean, or, equivalently, between the data’s dispersion and their mean.

# Verification of the dispersion is an aspect of data quality control: how well do the data accord to the expectations of our analytic approach?

# estimateDispersions: First, it estimates a dispersion value for each gene, then it fits a curve through the estimates. Finally, it assigns to each gene a dispersion value, using a choice between the per-gene estimate and the fitted value.

data = estimateDispersions( data, method="blind", sharingMode="fit-only" , fitType="parametric")

# Now, we can attempt to find differential expression:
res = nbinomTest( data, colnames( count_table )[1], colnames( count_table )[2] )

resp=subset(res,res$padj<0.1)

res_pvalue=resp[ order(resp$pval), ]

res_down=resp[ order( resp$foldChange, -resp$baseMean ), ]

res_up=resp[ order( -resp$foldChange, -resp$baseMean ), ]

write.csv(res_pvalue, file = sprintf("%sDESeq_%s_pvalue.csv", args[3], args[2]))

write.csv(res_up, file = sprintf("%sDESeq_%s_up.csv", args[3], args[2]))

write.csv(res_down, file = sprintf("%sDESeq_%s_down.csv", args[3], args[2]))

write.csv(res, file = sprintf("%sDESeq_%s_ALL.csv", args[3], args[2]))

