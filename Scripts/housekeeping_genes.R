## Load data
args <- commandArgs(TRUE)
#args <- c("1471-2229-10-160-s8.csv", "./")
count_table = read.csv(file = args[1], header=TRUE, row.names=1)

library( "DESeq" )


## Define data design
dataDesign = data.frame(
        row.names = colnames( count_table ),
        condition = colnames( count_table ),
        libType = rep("paired-end",ncol(count_table)))

conditions=dataDesign$condition

libraryType=dataDesign$libType

data = newCountDataSet(count_table, conditions)


## 1. NORMALIZATION
# estimate size factors (which is a representation of the differences in coverage between replicates)
# The effective library size information is called the size factors vector, since the package only needs to know the relative library sizesdata

data = estimateSizeFactors(data)

normalized=counts(data)
#normalized = counts(data, normalized=TRUE)

## Coeficient of variation
cv_calculation <- function(normalized){
      
	cv_table=c()
	rownames(cv_table)=c()
	
	for (i in 1:nrow(normalized)){
		#standart deviation
		st_dev=sd(normalized[i,])
		#mean
		u=mean(normalized[i,])
		#coeffecient of variation
		cv=st_dev/u
		
		cv_rowname=rownames(cv_table)
		cv_table=matrix(c(cv_table, cv), ncol=1, byrow=TRUE)
		rownames(cv_table)=c(cv_rowname, rownames(normalized)[i])

	}
	colnames(cv_table)=c("Coeficient_of_variation")
	return(cv_table)
}

cv_table=cv_calculation(normalized)
cv_table_up=cv_table[ order(Coeficient_of_variation), ]

write.csv(cv_table_up, file = sprintf("%housekeeping_genes_cv.and.deseq.csv", args[2]))
