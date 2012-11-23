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
entrezIds <- mget(genes_refseq, envir=org.Hs.egREFSEQ2EG)
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
### code chunk number 11: standardHyperGeo
###################################################
hgCutoff <- 0.01
params <- new("GOHyperGParams",
              geneIds=genesIds, #or gene_entrez (equal)
              universeGeneIds=entrezUniverse,
              annotation="org.Hs.eg.db",
              ontology="BP",
              pvalueCutoff=hgCutoff,
              conditional=FALSE,
              testDirection="over")

###################################################
### code chunk number 12: condHyperGeo
###################################################
paramsCond <- params

###################################################
### code chunk number 13: standardHGTEST
###################################################
hgOver <- hyperGTest(params)

###################################################
### code chunk number 15: summaryEx
###################################################
df <- summary(hgOver)

write.table(res_pvalue, file = sprintf("%sGOsummary_%s_pvalue.txt", args[3], arg[2]))
