## Load data
args <- commandArgs(TRUE)
count_table_read = read.table(file = args[1], header=TRUE, row.names=1)
count_table=count_table_read[1:(nrow(count_table_read)-5),]


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
# estimate the size factors (which is a representation of the differences in coverage between replicates)
# The effective library size information is called the size factors vector, since the package only needs to know the relative library sizesdata

data = estimateSizeFactors(data)

## Histogram of the ratios to double-check the size factors.
# calculate the gene-wise geometric means
#geomeans <- exp( rowMeans( log( counts(data) ) ) )

# choose a sample whose size factor estimate we want to check
#j <- 1

# The size factor was estimated as described above, so the following two lines give the same result.
#print( sizeFactors(data)[j] )
#print( median( ( counts(data)[,j]/geomeans )[ geomeans>0 ] )  )

# Plot a histogram of the ratios of which we have just taken the median:
#hist( log2( counts(data)[,j] / geomeans ), breaks=100 )

# This histogram should be unimodal, with a clear peak at the value of the size factor. Mark it in red:
#abline( v=log2( sizeFactors(cds)[ j ] ), col="red" )


## 2. VARIANCE

# The inference in DESeq relies on an estimation of the typical relationship between the data’s variance and their mean, or, equivalently, between the data’s dispersion and their mean.

# Verification of the dispersion is an aspect of data quality control: how well do the data accord to the expectations of our analytic approach?

# estimateDispersions: First, it estimates a dispersion value for each gene, then it fits a curve through the estimates. Finally, it assigns to each gene a dispersion value, using a choice between the per-gene estimate and the fitted value.

data = estimateDispersions( data, method="blind", sharingMode="fit-only" , fitType="local")

#This function calculates a variance stabilising transformations (VST) from the fitted dispersion-mean reltions and then transforms the count data (after normalization by division by the size factor), yielding a matrix of values which are now approximately homoskedastic. This is useful as input to statistical analyses requiring homoskedasticity.

vds=getVarianceStabilizedData(data)
library("RColorBrewer")
library("gplots")

## Heatmap of the genes with higher differential expression (first 30)
select = order(rowMeans(counts(data)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

## Heatmap after VST
png(filename= sprintf("%s%s.heatmap.with.vst.data.(30 genes).png", args[2], args[3]))
heatmap.2(vds[select,], col = hmcol, trace="none", margin=c(10, 6))
dev.off()

## Heatmap before VST
png(filename= sprintf("%s%s.heatmap.without.vst.data.(30 genes).png", args[2], args[3]))
heatmap.2(counts(data)[select,], col = hmcol, trace="none", margin=c(10,6))
dev.off()
