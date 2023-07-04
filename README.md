# DRIE
A drug repurposing method based on inhibition effect on gene regulatory network

# More about DRIE
The DRIE is a novel tool used to identify the potentail candidate drugs for cancers based on gene regulatory network.

# Getting started

## Step 1. Pre run the method for installation

You should ensure that you have the necessary system dependencies configured.

For Windows (8.1 / 10 / 11): Rtools should be installed to the system path.

The latest base R is recommended. The compatibility of the earlier version (v3.6.x) is under evaluation.
We use R version is [64-bit] d:\Program Files\R\R-3.6.0

## Step 2. Install the package
The dependency `EnrichmentBrowser`, `KEGGdzPathwaysGEO`, `KEGGandMetacoreDzPathwaysGEO` and `SPIA` are unavailable on the CRAN but available on [BioConductor](https://www.bioconductor.org/). So we need to install the BiocManager manually. 

``` r
if (!"BiocManager" %in% as.data.frame(installed.packages())$Package)
  install.packages("BiocManager")
BiocManager::install(c("EnrichmentBrowser", "KEGGdzPathwaysGEO","KEGGandMetacoreDzPathwaysGEO","SPIA"))
```

Then you can install the development version of DRIE from [GitHub](https://github.com/) with:

``` r
if (!"devtools" %in% as.data.frame(installed.packages())$Package)
  install.packages("devtools")
devtools::install_github("eshinesimida/DRIE")

```
## Examples

Below is a basic example that shows how to obtain candidate drugs of colorecatal cancer:

``` r
#Load require package
library(EnrichmentBrowser)
library(KEGGdzPathwaysGEO)
library(KEGGandMetacoreDzPathwaysGEO)
library(SPIA)
library(DRIE)


#1.Load the data of lung cancer
#--------GSE19188 differnential analysis-- lung cancer
data("GSE19188")
exprs_all <- exprs(GSE19188)
# Add the gene symbol
all.eset <- probe2gene(GSE19188)

before.norm <- assay(all.eset)
# Gene normalization
all.eset <- normalize(all.eset, norm.method="quantile")
after.norm <- assay(all.eset)

exprs_all1 <- data.frame(after.norm)
table(colData(all.eset)$Group)
colData(all.eset)$GROUP <- ifelse(colData(all.eset)$Group == "d", 1, 0)
normal <- length(which(colData(all.eset)$GROUP == '0'))
tumor <- length(which(colData(all.eset)$GROUP == '1'))
# Get the differential expression genes in limmar package
all.eset <- deAna(all.eset, padj.method="BH")
all_de <- rowData(all.eset, use.names=TRUE)
all_de <- data.frame(all_de)
tg <- all_de[order(all_de$ADJ.PVAL, decreasing = FALSE),]

tg3 <- all_de[all_de$ADJ.PVAL<0.05,]

#2. get disease-related genes of lung cancer from kegg pathways
library(DrugDiseaseNet)
library(graph)
gg<-keggGlobalGraph()
kegg_genes <- process(nodes(gg))
tg_3 <- tg3[intersect(rownames(tg3),kegg_genes),]
d <- tg_3$FC
names(d) <- rownames(tg_3)
d1 <- ifelse(d > 0,1,-1)
disease_genes <- d1
##------------------------------

#3.get lung cancer-related pathways from KEGG database
library(DrugDiseaseNet)
library(graph)
library(graphite)
library(Rgraphviz)
library(SPIA)
library(XML)


kgml.path=system.file("extdata/keggxml/hsa",package="SPIA")
mdir = kgml.path
paths <- dir(mdir,pattern=".xml")
genes_CRC <- c()
for(i in 1:8){
  gg <-try(parseKGML(paste(mdir,paths[i],sep="/")),TRUE)

  mapkG3<-KEGGpathway2Graph(gg,expandGenes=T)
  genes_CRC <- c(genes_CRC, nodes(mapkG3))


}
genes_CRC1 <- unique(genes_CRC)

gg<-keggGlobalGraph()
kegg_genes <- process(nodes(gg))
genes1 <- c('hsa:4790', 'hsa:5599', 'hsa:1839', 'hsa:253959', 'hsa:57148', 'hsa:57186', 'hsa:53373',
            'hsa:4486', 'hsa:844', 'hsa:845', 'hsa:10345', 'hsa:444', 'hsa:3270', 'hsa:255231', 'hsa:55283',
            'hsa:57192', 'hsa:2668', 'hsa:4485', 'hsa:10393', 'hsa:119504', 'hsa:246184', 'hsa:25847', 'hsa:25906',
            'hsa:29882', 'hsa:29945', 'hsa:51433', 'hsa:51434', 'hsa:51529', 'hsa:64682', 'hsa:8697', 'hsa:8881',
            'hsa:996', 'hsa:4088', 'hsa:10572', 'hsa:84883')

m1 <- setdiff(genes_CRC1 , genes1)
gg1 <- subGraph(m1,gg)
mapkG3<-gg1

library(igraph)
#g4 <- graph_from_graphnel(mapkG3)
adjmatrix<-as(mapkG3,"matrix")
globalGraphNonNeg<-mapkG3
mat<-adjmatrixNeg<-which(adjmatrix == -1,arr.ind=TRUE)
mat[,1]<-rownames(adjmatrix)[adjmatrixNeg[,1]]
mat[,2]<-rownames(adjmatrix)[adjmatrixNeg[,2]]

edgeData(globalGraphNonNeg,from=as.character(mat[,1]),
         to=as.character(mat[,2]),attr="weight")<-2
g4<-igraph.from.graphNEL(globalGraphNonNeg)

##------------------------------



