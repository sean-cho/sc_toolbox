tcgaGetNames = function(x){
  if(inherits(x,c("data.frame","matrix"))){stop("Data class cannot be matrix or data.frame.")}
  if(inherits(x,"list")){x = unlist(x)}
  if(!inherits(x,"character")){x = as.character(x)}
  if(nchar(x[1])>=8){
    return(substr(x,9,12))
  }
  else{
    return(substr(x,1,4))
  }
}

tcgaGetType = function(x){
  if(inherits(x,c("data.frame","matrix"))){stop("Data class cannot be matrix or data.frame.")}
  if(inherits(x,"list")){x = unlist(x)}
  if(!inherits(x,"character")){x = as.character(x)}
  if(nchar(x[1])>=8){
    return(substr(x,14,15))
  }
  else{
    return(substr(x,6,7))
  }
}

tcgaGetBoth  = function(x){
  if(inherits(x,c("data.frame","matrix"))){stop("Data class cannot be matrix or data.frame.")}
  if(inherits(x,"list")){x = unlist(x)}
  if(!inherits(x,"character")){x = as.character(x)}
  if(nchar(x[1])>=8){
    return(substr(x,9,15))
  }
  else{
    print("Vector appears to be in the right length")
    return(x)
  }
}

ggn = function(x){strsplit(x,"|",fixed=TRUE)[[1]][1]}

ggn.number = function(x){strsplit(x,"|",fixed=TRUE)[[1]][2]}

processTTB = function(x){
  dummy = rownames(x)
  rownames(x) = sapply(dummy,ggn.number)
  x$Gene = sapply(dummy,ggn)
  return(x)
}

getSigGene = function(x,dir="all",p=0.05,fc=1){
  if(dir=="all"){
    return(x$Gene[x$adj.P.Val<=p & abs(x$logFC)>=fc])
  }
  if(dir=="up"){
    return(x$Gene[x$adj.P.Val<=p & x$logFC>=fc])
  }
  if(dir=="down"){
    return(x$Gene[x$adj.P.Val<=p & x$logFC<=-fc])
  }
}
# 
# 
# tcgaEAnnot = sapply(rownames(tcgaE.TN),ggn.number)
# names(tcgaEAnnot) = sapply(rownames(tcgaE.TN),ggn)