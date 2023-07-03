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
library(DRIE)


#Load the data

#--------1.GSE19188 differnential analysis-- lung cancer
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

