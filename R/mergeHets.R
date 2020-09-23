##' Merge together heterozygous SVs that are very similar. This should
##' be done for SVs of the same type and genotype. When several SV pairs
##' can be merged in a cluster, only the best pair (best overlap) is
##' merged, hence this could be run several times (e.g. until the
##' size of the output is similar to the input's).
##' @title Merge heterozygous variants
##' @param svs GRanges with SV information
##' @param min.rol minimum reciprocal overlap to match variants.
##' @param max.ins.dist maximum distance for insertions to be clustered.
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @return an upated GRanges.
##' @author Jean Monlong
##' @keywords internal
mergeHets <- function(svs, min.rol=.9, max.ins.dist=1, ins.seq.comp=FALSE){
  if(length(svs)==0){
    return(svs)
  }
  if(length(unique(svs$type))>1){
    stop('mergeHets should be run separately for each type.')
  }

  ## overlap SVs
  ol.gr = prepareOl(svs, svs, min.rol=min.rol, max.ins.dist=max.ins.dist, ins.seq.comp=ins.seq.comp)
  ol.gr = annotateOl(ol.gr, min.ol=min.rol, method='reciprocal')
  
  ## If nothing, return input variants
  if(length(ol.gr)==0){
    logging::loginfo('No variants to merge')
    return(svs)
  }
  logging::loginfo(paste(length(ol.gr), 'pairs of variants to merge'))
  ## Select pairs to merge, best overlap first
  ol.gr = ol.gr[order(-ol.gr$olScore),]
  dup = duplicated(as.vector(rbind(ol.gr$queryHits, ol.gr$subjectHits)))
  dup = matrix(dup, 2)
  dup = colSums(dup)
  ol.gr = ol.gr[which(dup==0),]
  ## merge pairs
  s1 = GenomicRanges::start(svs[ol.gr$queryHits])
  s2 = GenomicRanges::start(svs[ol.gr$subjectHits])
  e1 = GenomicRanges::end(svs[ol.gr$queryHits])
  e2 = GenomicRanges::end(svs[ol.gr$subjectHits])
  starts = (s1 + s2)/2
  ends = (e1 + e2)/2
  chrs = as.character(GenomicRanges::seqnames(svs))[ol.gr$queryHits]
  svs.merged = GenomicRanges::GRanges(chrs, IRanges::IRanges(starts, ends))
  ## Merge columns
  for(coln in colnames(GenomicRanges::mcols(svs))){
    if(coln == 'ac'){
      svs.merged$ac = svs$ac[ol.gr$queryHits] + svs$ac[ol.gr$subjectHits]
    } else if(coln == 'qual'){ # Should we use the average quality ???
       svs.merged$qual = (svs$qual[ol.gr$queryHits] + svs$qual[ol.gr$subjectHits])/2
    } else if(coln == 'ref'){
      svs.merged$ref = svs$ref[ol.gr$queryHits]
    } else if(coln == 'alt'){
      svs.merged$alt = svs$alt[ol.gr$queryHits]
    } else if(coln == 'type'){
      svs.merged$type = svs$type[ol.gr$queryHits]
    } else if(coln == 'size'){
      svs.merged$size = (svs$size[ol.gr$queryHits] + svs$size[ol.gr$subjectHits]) / 2
    }
  }
  ## Remove merged pairs and add new SVs
  svs = svs[setdiff(1:length(svs), unique(c(ol.gr$queryHits, ol.gr$subjectHits)))]
  svs = c(svs, svs.merged)
  return(svs)
}
