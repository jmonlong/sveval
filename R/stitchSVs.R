##' Stitch together SVs that are close together. This should be done for SVs of the same
##' type and genotype. When several SV pairs can be stitched in a cluster, only the best pair (closest to each other) is merged, hence this should be run several times (e.g. until the
##' size of the output is similar to the input's).
##' @title Stitch fragmented SVs
##' @param svs GRanges with SV information
##' @param stitch.dist the maximum distance allowed for two SVs to be stitched.
##' @return an upated GRanges.
##' @author Jean Monlong
##' @keywords internal
stitchSVs <- function(svs, stitch.dist=20){
  if(length(svs)==0){
    return(svs)
  }
  if(length(unique(svs$type))>1){
    stop('stitchSVs should be run separately for each type.')
  }
  ## Overlap against itselfGenomicRanges::
  ol = GenomicRanges::findOverlaps(svs, svs, maxgap=stitch.dist)
  ol = as.data.frame(ol)
  ## Remove identity and redundant pairs
  ol = ol[which(ol$queryHits < ol$subjectHits),]
  ## Compute distance
  ol$d = GenomicRanges::distance(svs[ol$queryHits], svs[ol$subjectHits])
  ## fix distance=0 for immediately adjacent regions
  ol.any = IRanges::overlapsAny(svs[ol$queryHits], svs[ol$subjectHits])
  ol$d = ifelse(ol$d==0 & !ol.any, 1, ol$d)
  ## Filter overlapping pairs, we want to stitch fragmented calls
  ## and we assume that fragmented pieces don't overlap
  ol = ol[which(ol$d>0),]
  ## If nothing, return input variants
  if(nrow(ol)==0){
    logging::loginfo('no variants to stitch')
    return(svs)
  }
  logging::loginfo(paste(nrow(ol), 'pairs of variants to stitch'))
  ## Select pairs to stitch, nearest to each other first
  ol = ol[order(ol$d),]
  dup = duplicated(as.vector(rbind(ol$queryHits, ol$subjectHits)))
  dup = matrix(dup, 2)
  dup = colSums(dup)
  ol = ol[which(dup==0),]
  ## Stitch pairs
  s1 = GenomicRanges::start(svs[ol$queryHits])
  s2 = GenomicRanges::start(svs[ol$subjectHits])
  e1 = GenomicRanges::end(svs[ol$queryHits])
  e2 = GenomicRanges::end(svs[ol$subjectHits])
  starts = ifelse(s1<s2, s1, s2)
  ends = ifelse(e1>s2, e1, e2)
  chrs = as.character(GenomicRanges::seqnames(svs))[ol$queryHits]
  svs.stitch = GenomicRanges::GRanges(chrs, IRanges::IRanges(starts, ends))
  ## Merge columns
  for(coln in colnames(GenomicRanges::mcols(svs))){
    if(coln == 'ac'){
      svs.stitch$ac = svs$ac[ol$queryHits]
    } else if(coln == 'qual'){ # Should we use the average quality ???
       svs.stitch$qual = (svs$qual[ol$queryHits] + svs$qual[ol$subjectHits])/2
    } else if(coln == 'alt'){
      svs.stitch$alt = paste(svs$alt[ol$queryHits], svs$alt[ol$subjectHits], sep='NNN')
    } else if(coln == 'type'){
      svs.stitch$type = svs$type[ol$queryHits]
    } else if(coln == 'size'){
      svs.stitch$size = svs$size[ol$queryHits] + svs$size[ol$subjectHits]
    } else if(coln == 'ref.cov'){
      svs.stitch$ref.cov = svs$ref.cov[ol$queryHits] + svs$ref.cov[ol$subjectHits]
    } else if(coln == 'alt.cov'){
      svs.stitch$alt.cov = svs$alt.cov[ol$queryHits] + svs$alt.cov[ol$subjectHits]
    }
  }
  ## Remove stitched pairs and add new SVs
  svs = svs[setdiff(1:length(svs), unique(c(ol$queryHits, ol$subjectHits)))]
  svs = c(svs, svs.stitch)
  return(svs)
}
