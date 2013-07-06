## Load count data
args <- commandArgs(TRUE)
library( "DESeq" )

args<-c("LM1_HeLa_20_control.txt",".","./")
count_table_read_LM1 = read.table(file = args[1], header=TRUE, row.names=1)
#count_table_read_LM2 = read.table(file = args[2], header=TRUE, row.names=1)

library( "DESeq" )


## Define data design
# LM1
dataDesign_LM1 = data.frame(
        row.names = colnames( count_table_read_LM1 ),
        condition = colnames( count_table_read_LM1 ),
        libType = rep("paired-end",ncol(count_table_read_LM1)))
conditions_LM1=dataDesign_LM1$condition
libraryType_LM1=dataDesign_LM1$libType

n_LM1=nrow(count_table_read_LM1)
count_table_LM1=count_table_read_LM1[-((n_LM1-4):n_LM1),]
data = newCountDataSet(count_table_LM1, conditions_LM1)

# LM2
dataDesign_LM2 = data.frame(
        row.names = colnames( count_table_read_LM2 ),
        condition = colnames( count_table_read_LM2 ),
        libType = rep("paired-end",ncol(count_table_read_LM2)))
conditions_LM2=dataDesign$condition
libraryType_LM2=dataDesign$libType

n_LM2=nrow(count_table_read_LM2)
count_table_LM2=count_table_read_LM2[-((n_LM2-4):n_LM2),]
data_LM2 = newCountDataSet(count_table_LM2, conditions_LM2)


## Load DEG
all_pvalue_files=list.files(args[3], "*pvalue.csv")
files_LM1=all_pvalue_files[grep("LM1", all_pvalue_files)]
files_LM2=all_pvalue_files[grep("LM2", all_pvalue_files)]




## Get counts for DEG 
get_genes_per_cond = function(files_DEG, count_tab){
	deg_genes=c()
	for( i in 1:length(files_DEG)){
		table = read.csv(file=sprintf("%s%s", args[3], files_DEG[i]), head=TRUE, sep=",", row.names=1)
		genes_col=table$id
		for( j in 1:length(genes_col)){
			pos_gene = genes_col[j] == row.names(count_tab)
			line_count_tab = as.matrix(count_tab[pos_gene,])
		}
		deg_genes= rbind(deg_genes,line_count_tab)
	}
}

files_DEG=files_LM1
count_tab=count_table_LM1
get_genes_per_cond(files_LM1, count_table_LM1)
