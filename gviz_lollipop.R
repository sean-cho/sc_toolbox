# Readme ------------------------------------------------------------------

### Program: Gviz_lollipop.R

### Description
# This code adds a function to make (pseudo)-lollipop plots from Illumina HM450K 
# data to the Gviz package from bioconductor.
# The lollipop fills, by default, on a 0 - 1 Beta value scale
# Dependencies: Gviz, minfi, GenomicRanges, BSGenome, scales
#               FDb.InfiniumMethylation.hg19
# Genome used is hg19

### IMPORTANT
# 1. If running from shell, remove examples section before running R CMD BATCH
# 2. Remember to execute code from Exit section to reset options() to default

# Arguments ---------------------------------------------------------------

### USER DEFINE:
# load('path_to/your_data.rda')

### Input data
# rgset: Should be a MethylSet/RGset with appropriate preprocessing
rgset = tn_fn
# betaset:
# If betaset is NULL, will obtain using getBeta on rgset
# You may specify your own, if you have prefiltered for probes
betaset = NULL

### Required arguments

# poi: {character} probe of interest (only one is allowed)
poi = 'cg05558259'

### Predefined arguments

# adj_bases: {integer} How many bases around the probe to extend
adj_bases = 500L
# y_padding: {numeric} The amount of padding on the y-axis on the Gviz track
y_padding = 0.05
# pointsize: {numeric} Defines gpar cex for size of lollipops
pointsize = 0.5
# track_name: {character} Name of the track
track_name = 'CpG Lollipop'

# Setup -------------------------------------------------------------------

require(scales)
require(BSgenome.Hsapiens.UCSC.hg19)
require(minfi)
require(GenomicRanges)
require(Gviz)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)

# Stop the code if an error is seen
default_error <- getOption('error')
options(error = utils::dump.frames)

# Define objects (run once) -----------------------------------------------

# Get 450k annotation information
hm450 <- getAnnotation(rgset)
if(is.null(betaset)) betaset <- getBeta(rgset)

hm450 <- hm450[rownames(betaset),]
hm450$probe_id <- rownames(hm450)
dummy_seqinfo <- seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene)[paste0('chr',c(1:22,'X','Y'))]

# Convert data.frame to GenomicRanges
hr450 <- with(hm450, GRanges(seqnames = chr, 
                             ranges = IRanges(start = pos, width = 1), strand = strand,
                             probe_id = probe_id,
                             seqinfo = dummy_seqinfo))
names(hr450) <- hr450$probe_id
rm(dummy_seqinfo)

# Sub-workflows -----------------------------------------------------------

# annotate_cg finds all nearby CpGs and probes for downstream functionality
annotate_cg <- function(poi, probe_region = adj_bases, plot_region = 1L){
  require(BSgenome.Hsapiens.UCSC.hg19)
  hr450 <- with(hm450, GRanges(seqnames = chr, 
                               IRanges(start = pos, width = 1), 
                               strand = strand,
                               probe_id = Name))
  
  probehm <- hm450[poi,]
  probehr <- subset(hr450, probe_id == poi)
  probe_flank <- flank(probehr, width = probe_region, both = TRUE)
  nearby_probes <- subsetByOverlaps(hr450, probe_flank, type = 'any')
  cat('Found', length(nearby_probes), 'adjacent probes.\n')
  
  strandedness <- as.character(strand(probehr))
  
  hgenome <- BSgenome.Hsapiens.UCSC.hg19
  region_loc <- range(probe_flank)
  region_seq <- getSeq(hgenome, region_loc)
  if(strandedness == '-') region_seq <- reverseComplement(region_seq)
  region_cg <- vmatchPattern('CG', region_seq)
  cat('Found', length(region_cg[[1]]), 'adjacent CpG dinucleotides.\n')
  
  cg_starts = start(region_loc) + start(region_cg[[1]]) - 1
  cg_ranges = GRanges(seqnames = seqnames(region_loc), 
                      ranges = IRanges(start = cg_starts, width = 1),
                      strand = strand(region_loc))
  cg_ranges$Illumina = FALSE
  cg_ranges[findOverlaps(cg_ranges, nearby_probes, type = 'any')@queryHits]$Illumina = TRUE
  
  cg_ranges$probe_id = NA_character_
  cg_ranges[findOverlaps(cg_ranges, nearby_probes, type = 'any')@queryHits]$probe_id =
    nearby_probes[findOverlaps(nearby_probes, cg_ranges, type = 'any')@queryHits]$probe_id
  
  output <- list()
  output$probes <- nearby_probes
  output$cpgs <- cg_ranges
  return(output)
}

# lollicol scales the color of beta-values with 0 = white, 1 = black
lollicol <- function(x){
  y <- 1 - x
  return(rgb(red = y, green = y,blue = y))
}

# Workhorse ---------------------------------------------------------------

cpgs <- annotate_cg(poi = poi, probe_region = adj_bases)

ct_input <- t(betaset[cpgs$probes$probe_id,,drop = FALSE])
lower_bound <- 0 + y_padding
upper_bound <- 1 - y_padding

ct = CustomTrack(plottingFunction = function(GdObject, prepare) {
  if(!prepare) {
    yscales <- scales::rescale(x = 1:nrow(ct_input), 
                               to = c(lower_bound, upper_bound))
    for(i in 1:nrow(ct_input)){
      grid.points(x = start(cpgs$probes), y = rep(yscales[i], ncol(ct_input)), pch = 16,
                  gp = gpar(col = lollicol(ct_input[i,]), 
                            cex = pointsize))
      grid.points(x = start(cpgs$probes), y = rep(yscales[i], ncol(ct_input)),
                  gp = gpar(cex = pointsize, 
                            col = 'black'))
    }
  }
  return(invisible(GdObject))}, name = track_name)

# Usage -------------------------------------------------------------------

coi <- as.character(seqnames(hr450[poi,]))
roi <- range(cpgs$cpgs)
plot_from <- start(roi)
plot_to <- end(roi)

# Without annotation
plotTracks(ct, chromosome = coi, from = plot_from, to = plot_to)

# With annotation
# Make annotation track
gviz_cg <- AnnotationTrack(range = cpgs$cpgs, shape = 'box', name = 'CpGs')

plotTracks(list(gviz_cg, ct),
           chromosome = coi, from = plot_from, to = plot_to,
           sizes = c(1, 20))

# Exit --------------------------------------------------------------------

options(error = default_error)

# Example(s) --------------------------------------------------------------

gviz_id <- IdeogramTrack(genome = 'hg19')
gviz_ax <- GenomeAxisTrack()

gviz_set <- list(gviz_id, gviz_ax, gviz_cg, ct)
plotTracks(gviz_set,
           chromosome = coi, from = plot_from, to = plot_to,
           sizes = c(1, 2, 1, 20))