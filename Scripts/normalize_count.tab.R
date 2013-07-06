args <- commandArgs(TRUE)
library( "DESeq" )

count_table = read.table(file = args[1], header=TRUE, row.names=1)

## Define data design
dataDesign = data.frame(
        row.names = colnames( count_table ),
        condition = colnames( count_table ),
        libType = rep("paired-end",ncol(count_table)))
conditions=dataDesign$condition
libraryType=dataDesign$libType

data = newCountDataSet(count_table, conditions)
data = estimateSizeFactors(data)
normalized = counts(data, normalized=TRUE)
write.table(normalized, file=sprintf("%s_normalized", args[1]))
