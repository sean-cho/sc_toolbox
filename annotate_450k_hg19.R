annotate_450k_hg19 <- function(gene = TRUE, tss = TRUE, tx = TRUE, 
                          classify_type = TRUE, verbose = TRUE){
  
  require(FDb.InfiniumMethylation.hg19)
  
  # Get genomic coordinates for 450k probes
  meta_450k <- get450k()
  
  # Convert to data.frame output
  metadf_450k <- as.data.frame(meta_450k)
  
  # Map to genes
  if(gene & !tx) {
    if(verbose) cat('Getting gene information.\n')
    gene_450k <- getNearestGene(meta_450k)
    if(!all(rownames(metadf_450k) == rownames(gene_450k))) stop('Metadata and gene probeids do not match.')
    metadf_450k$gene <- gene_450k$nearestGeneSymbol
  }
  if(tx) {
    if(verbose) cat('Getting gene & tx information.\n')
    tx_450k <- getNearestTranscript(meta_450k)
    if(!all(rownames(metadf_450k) == rownames(tx_450k))) stop('Metadata and tx probeids do not match.')
    metadf_450k$gene <- tx_450k$nearestGeneSymbol
    metadf_450k$tx <- tx_450k$nearestTranscript
  }
  if(tss) {
    if(verbose) cat('Getting tss information.\n')
    tss_450k <- getNearestTSS(meta_450k)
    if(!all(rownames(metadf_450k) == rownames(tss_450k))) stop('Metadata and tss probeids do not match.')
    metadf_450k$tss_gene <- tss_450k$nearestGeneSymbol
    metadf_450k$tss_tx <- tss_450k$nearestTranscript
  } 
  
  # Classify probe type
  if(classify_type) {
    cat('Classifying probe types.\n')
    metadf_450k$type <- 'CpG'
    metadf_450k$type[grep('rs', rownames(metadf_450k))] <- 'Control'
    metadf_450k$type[grep('ch', rownames(metadf_450k))] <- 'non-CpG'
  }
  
  return(metadf_450k)
}

# Example
# ann_450k <- annotate_450k()
