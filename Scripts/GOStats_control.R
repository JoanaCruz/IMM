args <- commandArgs(TRUE)

###################################################
### Setup
###################################################
library("org.Hs.eg.db")
library("GO.db")
library("annotate")
library("genefilter")
library("GOstats")
library("RColorBrewer")
library("xtable")


###################################################
### Load data
###################################################
data <- read.csv(file=args[1], head=TRUE, sep=",", row.names=1)
file_name = strsplit(args[1],'_')

if(length(data$id)==0){
	write("No differencially expressed genes", file = sprintf("%sNO_GOsummary_%s_%s.txt", args[2], file_name[[1]][2], file_name[[1]][4]))
} else {

genes_refseq=as.character(data$id)


###################################################
### Filtering-noEntrez
###################################################
## Remove genes that have no entrezGene id
entrezIds <- mget(genes_refseq, envir=org.Hs.egREFSEQ2EG, ifnotfound=NA)
genes_refseq <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]


###################################################
### Filtering-noGO
###################################################
## Convert data to EntrezId
genes_entrez=c()
for(i in 1:length(genes_refseq)){
	genes_entrez=c(genes_entrez,get(genes_refseq[i],org.Hs.egREFSEQ2EG))
}

## Remove genes with no GO mapping
haveGo <- sapply(mget(genes_entrez, org.Hs.egGO),
                 function(x) {
                     if (length(x) == 1 && is.na(x)) 
                       FALSE 
                     else TRUE
                 })
numNoGO <- sum(!haveGo)
genes_entrez <- genes_entrez[haveGo]


###################################################
### defineGeneUniverse
###################################################
## Convert data to Refseq
genes_refseq=c()
for(i in 1:length(genes_entrez)){
	genes_refseq=c(genes_refseq,get(genes_entrez[i],org.Hs.egREFSEQ))
}

## Relation between Refseq and Entrez Id
genesIds <- mget(genes_refseq, org.Hs.egREFSEQ2EG)

## Define gene universe (relation between all Refseq and Entrez Id)
entrezUniverse <- as.list(org.Hs.egREFSEQ2EG)


###################################################
### HyperGeo
###################################################
hgCutoff <- 0.05

ontologies = c('MF', 'BP', 'CC')
for (ont in ontologies) {
	params <- new("GOHyperGParams",
             	geneIds=genesIds, #or gene_entrez (equal)
              	universeGeneIds=entrezUniverse,
              	annotation="org.Hs.eg.db",
              	ontology=ont,
              	pvalueCutoff=hgCutoff,
              	conditional=FALSE,
              	testDirection="over")

paramsCond <- params

###################################################
### HGTEST
###################################################
hgOver <- hyperGTest(params)


###################################################
### summary
###################################################
df <- summary(hgOver)

goterms <- df[,1]
i=1
genes_ids <- c()

while (i <= length(goterms)){
	genes_ids <- c(genes_ids ,geneIdsByCategory(hgOver, goterms[i]))
	i <- i+1
}


## Convert data to EntrezId
for(i in 1:length(genes_ids)){
	genes_refseq=c()
	if(length(genes_ids[[i]]>0)){
		for(j in 1:length(genes_ids[[i]])){
			all_refseq = get(genes_ids[[i]][j],org.Hs.egREFSEQ)
			mrna_refseq_table = grep("NM", as.character(all_refseq))
			mrna_refseq=all_refseq[mrna_refseq_table]
			genes_refseq = c(genes_refseq, mrna_refseq)
		}
	}
	genes_ids[[i]]=genes_refseq
}

## Add column with genes (RefSeq annotation) 
df$genes_RefSeq <- genes_ids

df$genes_RefSeq <- sapply(df$genes_RefSeq, FUN = paste, collapse = ", ")
print(strsplit(args[1],'_'))

## Write table
write.csv(df, file = sprintf("%sGOsummary_%s_%s_%s.csv", args[2], file_name[[1]][2], ont, file_name[[1]][4]))

}
}