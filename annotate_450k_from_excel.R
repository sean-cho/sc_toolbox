annotate_450k_from_excel <- function(data_matrix) {
  require(minfi)
  require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  
  MS <- new('MethylSet', Meth = data_matrix, Unmeth = data_matrix,
            annotation = c('array' = 'IlluminaHumanMethylation450k',
                           'annotation' = 'ilmn12.hg19'))
  ann <- getAnnotation(MS)
  return(ann)
}

### Example
# mset = read.xls('path/to/data.xls', stringsAsFactors = FALSE)
# rownames(mset) = mset$TargetID
# keep = grep('AVG_Beta', colnames(mset), value = TRUE)
# betaset = subset(mset, select = keep)
# ann_450k = annotate_450k_from_excel(betaset)

# It would also be acceptable to use a pre-filtered betaset, after filtering against
# probes with poor call rates