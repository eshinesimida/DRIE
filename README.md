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
