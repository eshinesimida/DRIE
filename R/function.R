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