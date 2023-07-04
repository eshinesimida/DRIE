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
source('function.R')


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
mapkG3 <- get_pathway()

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

#-----4.get drug-induced genes from CMap database
drug_genes <- read.csv('drug_gene_all1.csv', header = T, stringsAsFactors = F)
drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'MCF7'),]
#drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'HL60'),]
#drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'PC3'),]
drug_genes_MCF7_1 <- drug_genes_MCF7[drug_genes_MCF7$fdr<0.05,]
MCF1 <- as.character(unique(drug_genes_MCF7_1$drug))
##------------------------------

#-----5. run the DRIE method based on lung cancer
drugs <- c()
scores1 <- c()
p_name <- c()
pb <- txtProgressBar(min = 0, max = length(MCF1), style = 3)

for(ii in 1:length(MCF1)){
  setTxtProgressBar(pb, ii)
  i1 <- MCF1[ii]
  cat('i=',i1,'\n')
  drug_gene <- drug_genes_MCF7_1[which(drug_genes_MCF7_1$drug == i1),]
  fc <- drug_gene$fc
  names(fc) <- drug_gene$ENTREZID
  drug_gene1 <- paste('hsa:',drug_gene$ENTREZID, sep = '')

  a <- length(intersect(nodes(mapkG3), drug_gene1))
  d_gene <- paste('hsa:',names(disease_genes),sep='')
  d_gene1 <- intersect(nodes(mapkG3), d_gene)
  if(a > 0 && length(d_gene1) > 0){

    diff_gene <- d_gene1
    diff_drug <- intersect(nodes(mapkG3), drug_gene1)

    distance1 <- c()
    score <- c()
    n = 0
    for(i in diff_gene){
      n = n+1
      cat('n=',n,'\n')
      for(j in diff_drug){
        I <- get_r2(g4, j ,i)
        score <- c(score,I)

      }
    }
    B <- matrix(score, nrow = length(diff_drug), byrow = FALSE, dimnames = list(diff_drug,
                                                                                diff_gene))

    B1 <- apply(B,2,sum)
    B2 <- c()
    for(m in B1){
      if(m >0){
        value <- 1
      }else if(m==0){
        value <- 0
      }else{
        value <- -1
      }
      B2 <- c(B2, value)
    }
    names(B2) <- names(B1)
    #B2 <- ifelse(B1< 0,-1,1)
    names(B2) <- sapply(names(B2),process)
    drug_gene_all <- B2
    disease_gene_all <- disease_genes[names(drug_gene_all)]
    consistence <- data.frame(drug_gene = drug_gene_all, disease_gene = disease_gene_all)
    #scores <- length(which(apply(consistence, 1, sum) == 0))
    scores <- sum(-consistence$drug_gene*consistence$disease_gene)
  }else{
    scores <- 0

  }
  drugs <- c(drugs, i1)
  scores1 <- c(scores1, scores)
  #p_name <- c(p_name, p)
}

dataframe_drug <- data.frame(drug = MCF1, score = scores1, stringsAsFactors = F)
##------------------------------
##------------------------------



