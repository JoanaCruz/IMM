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
data <- read.table(args[1], header=T, sep=" ", row.names=1)
genes_refseq=as.character(data$id)

###################################################
### Filtering-noEntrez
###################################################
## Remove genes that have no entrezGene id
entrezIds <- mget(genes_refseq, envir=org.Hs.egENSEMBL2EG, ifnotfound=NA)
genes_refseq <- names(entrezIds)[sapply(entrezIds, function(x) !is.na(x))]


###################################################
### Filtering-noGO
###################################################
## Convert data to EntrezId
genes_entrez=c()
for(i in 1:length(genes_refseq)){
	genes_entrez=c(genes_entrez,get(genes_refseq[i],org.Hs.egENSEMBL2EG))
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
	genes_refseq=c(genes_refseq,get(genes_entrez[i],org.Hs.egENSEMBL))
}

## Relation between Refseq and Entrez Id
genesIds <- mget(genes_refseq, org.Hs.egENSEMBL2EG)

## Define gene universe (relation between all Refseq and Entrez Id)
entrezUniverse <- as.list(org.Hs.egENSEMBL2EG)


###################################################
### HyperGeo
###################################################
hgCutoff <- 0.01

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

write.csv(df, file = sprintf("%sGOsummary_%s_%s.txt", args[3], ont, args[2]))
}
