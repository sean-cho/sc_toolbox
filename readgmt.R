## Reads a GMT file from msigdb and returns
## a list of [1] gene set names and [2] gene sets

parsegmt <- function(x){
    n <- length(x)
    return(x[3:n])
    }


readgmt <- function(gmtfile){
  gmtlist <- scan(gmtfile, what = 'character', sep = '\n')
  gmtlist <- strsplit(gmtlist, '\t')
  gsname <- sapply(gmtlist, '[[', 1)
  gsgenes <- lapply(gmtlist, parsegmt)
  names(gsgenes) <- gsname
  out <- list(gsnames = gsname,
              genesets = gsgenes)
  return(out)
  }