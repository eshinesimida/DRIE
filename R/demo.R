#disease gene
library(EnrichmentBrowser)
library(KEGGdzPathwaysGEO)
library(KEGGandMetacoreDzPathwaysGEO)
library(SPIA)


#function

get_r2 <- function(g4,from,to){
  
  #pathway <- prcocess_g(mapkG3)
  if(from == to){
    FC <- fc[strsplit(from,':')[[1]][2]]
    #if(FC >0){
    # y = 1
    # }else{
    # y = -1
    #}
    y <- FC
    
  }else{
    
    A <- get.shortest.paths(g4, from ,to)
    
    if(length(A$vpath[[1]]) == 1){
      y = 0
    }else{
      edge1 <- mapkG3@edgeData@data#边的类型
      L <- length(A$vpath[[1]])-1#最短距离
      E_type <- c()
      nodes <- A$vpath[[1]]
      for(i in 1:L){
        
        e <- edge1[paste(names(nodes[[i]]),'|',names(nodes[[i+1]]),sep = '')]
        E_type <- c(E_type,e[[1]]$weight)
      }
      
      FC <- fc[strsplit(from,':')[[1]][2]]
      S <- cumprod(E_type)
      s <- tail(S, 1)
      I <- FC*s/(L+1)
      #if(I > 0){
      #  y <- 1
      #}else{
      #  y <- -1
      #}
      y <- I
      
    }
    
    
  }
  y
}

process <- function(x){
  y <- gsub('hsa:','',x)
  y
}
#----------------------------------------
#--------1.GSE6956 differnential analysis--

all_de <- read.csv('Breast cancer_degs.csv', header = T, stringsAsFactors = F)

all_de <- read.csv('BC_gse31448.csv', header = T, stringsAsFactors = F)


all_de <- read.csv('BC_GSE42568.csv', header = T, stringsAsFactors = F)

all_de <- read.csv('BC_GSE29044.csv', header = T, stringsAsFactors = F)




#tg3 <- all_de[all_de$P.Value<0.05,]

tg3 <- all_de


genes <- c()
FC <- c()
pvalue <- c()
for(i in 1:length(tg3$Gene.ID)){
  m <- strsplit(tg3[i,]$Gene.ID,'///')[[1]]
  FC1 <- rep(tg3[i,]$logFC,length(m))
  pvalue1 <- rep(tg3[i,]$adj.P.Val, length(m))
  
  genes <- c(genes, m)
  FC <- c(FC, FC1)
  pvalue <- c(pvalue, pvalue1)
}

tg_all <- data.frame(genes = genes, FC = FC, pvalue = pvalue)




#
library(DrugDiseaseNet)
library(graph)
gg<-keggGlobalGraph()

kegg_genes <- process(nodes(gg))



tg <- read.csv('COAD.csv', header = T, stringsAsFactors = F)

tg <- tg[!duplicated(tg$ENTREZID),]
tg5 <- tg[tg$adj.P.Val<0.05,]

rownames(tg5) <- tg5$ENTREZID


tg_all <- tg_all[!duplicated(tg_all$genes),]
rownames(tg_all) <- tg_all$genes

tg_3 <- tg_all[intersect(tg_all$genes,kegg_genes),]

d <- tg_3$FC

names(d) <- rownames(tg_3)
d1 <- ifelse(d > 0,1,-1)
disease_genes <- d1

##------------------------------
#------------------------------

#crc-related pathways

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
for(i in 1:9){
  gg <-try(parseKGML(paste(mdir,paths[i],sep="/")),TRUE)
  
  mapkG3<-KEGGpathway2Graph(gg,expandGenes=T)
  genes_CRC <- c(genes_CRC, nodes(mapkG3))
  
  
}
genes_CRC1 <- unique(genes_CRC)

#MCF1 <- unique(drug_genes_MCF7_1$drug)

drugs <- c()
scores1 <- c()
p_name <- c()

library(DrugDiseaseNet)
library(graph)
gg<-keggGlobalGraph()

kegg_genes <- process(nodes(gg))

genes1 <- c('hsa:4790', 'hsa:5599', 'hsa:10393', 'hsa:119504', 'hsa:246184', 'hsa:25847', 
            'hsa:25906', 'hsa:29882', 'hsa:29945', 'hsa:51433', 'hsa:51434', 'hsa:51529', 
            'hsa:64682', 'hsa:8697', 'hsa:8881', 'hsa:996', 'hsa:4088', 'hsa:10572', 
            'hsa:84883', 'hsa:29964', 'hsa:4007', 'hsa:102723796', 'hsa:7088', 'hsa:7089',
            'hsa:7090', 'hsa:7091', 'hsa:79816', 'hsa:1499', 'hsa:2487', 'hsa:8840', 
            'hsa:4919', 'hsa:4920', 'hsa:54894', 'hsa:84133', 'hsa:55366', 'hsa:59352',
            'hsa:8549', 'hsa:284654', 'hsa:340419', 'hsa:343637', 'hsa:84870', 'hsa:25776',
            'hsa:1501', 'hsa:147495', 'hsa:164284', 'hsa:440193', 'hsa:342371', 'hsa:6310', 
            'hsa:1839', 'hsa:9166', 'hsa:100653049', 'hsa:125115', 'hsa:147183', 'hsa:162605', 
            'hsa:192666', 'hsa:25984', 'hsa:342574', 'hsa:353288', 'hsa:3857', 'hsa:3859', 
            'hsa:3860', 'hsa:3861', 'hsa:3866', 'hsa:3868','hsa:3872', 'hsa:3880', 'hsa:3881', 'hsa:3882', 'hsa:3883', 'hsa:3884', 'hsa:3885', 'hsa:3886', 'hsa:390792', 'hsa:54474', 'hsa:8687', 'hsa:8688', 'hsa:8689', 'hsa:7031')

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


#--------------------------------
#---------------------------------------
#----------------CMap dataset-----------------

#-------------------------------------------------
#-----9. MCF7, FDR<0.01-------------------------------------

drug_genes <- read.csv('drug_gene_all1.csv', header = T, stringsAsFactors = F)

drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'MCF7'),]


drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'HL60'),]
drug_genes_MCF7 <- drug_genes[which(drug_genes$cell_line == 'PC3'),]

drug_genes_MCF7_1 <- drug_genes_MCF7[drug_genes_MCF7$fdr<0.05,]



MCF1 <- as.character(unique(drug_genes_MCF7_1$drug))




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



