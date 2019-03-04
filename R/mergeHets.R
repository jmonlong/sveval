##' Merge together heterozygous SVs that are very similar. This should
##' be done for SVs of the same type and genotype. When several SV pairs
##' can be merged in a cluster, only the best pair (best overlap) is
##' merged, hence this should be run several times (e.g. until the
##' size of the output is similar to the input's).
##' @title Merge heterozygous variants
##' @param svs GRanges with SV information
##' @param min.rol minimum reciprocal overlap to match variants.
##' @param max.ins.gap maximum distance for insertions to be clustered.
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @return an upated GRanges.
##' @author Jean Monlong
##' @keywords internal
mergeHets <- function(svs, min.rol=.9, max.ins.gap=1, ins.seq.comp=FALSE){
  if(length(unique(svs$type))>1){
    stop('mergeHets should be run separately for each type.')
  }
  if(svs$type[1] == 'INS'){
    ol = olInsertions(svs, svs, max.ins.gap=max.ins.gap, ins.seq.comp=ins.seq.comp)
    ol = ol$ol
    ol = ol[which(ol$queryHits < ol$subjectHits),]
    rol.call = svs$size[ol$call.idx] / ol$call.cov
    rol.truth = svs$size[ol$truth.idx] / ol$truth.cov
    ol$rol = ifelse(rol.call<rol.truth, rol.call, rol.truth)
  } else {
    ol = GenomicRanges::findOverlaps(svs, svs)
    ol = as.data.frame(ol)
    ## Remove identity and redundant pairs
    ol = ol[which(ol$queryHits < ol$subjectHits),]
    ## Compute reciprocal overlap
    ol$qsw = GenomicRanges::width(GenomicRanges::pintersect(svs[ol$queryHits], svs[ol$subjectHits]))
    ol$qw = GenomicRanges::width(svs[ol$queryHits])
    ol$sw = GenomicRanges::width(svs[ol$subjectHits])
    ol$rol = ol$qsw / ifelse(ol$qw > ol$sw, ol$qw, ol$sw)
  }
  ## Filter overlapping pairs with low reciprocal overlap
  ol = ol[which(ol$rol>min.rol),]
  ## If nothing, return input variants
  if(nrow(ol)==0){
    return(svs)
    }
  ## Select pairs to merge, best overlap first
  ol = ol[order(-ol$rol),]
  dup = duplicated(as.vector(rbind(ol$queryHits, ol$subjectHits)))
  dup = matrix(dup, 2)
  dup = colSums(dup)
  ol = ol[which(dup==0),]
  ## merge pairs
  s1 = GenomicRanges::start(svs[ol$queryHits])
  s2 = GenomicRanges::start(svs[ol$subjectHits])
  e1 = GenomicRanges::end(svs[ol$queryHits])
  e2 = GenomicRanges::end(svs[ol$subjectHits])
  starts = (s1 + s2)/2
  ends = (e1 + e2)/2
  chrs = as.character(GenomicRanges::seqnames(svs))[ol$queryHits]
  svs.merged = GenomicRanges::GRanges(chrs, IRanges::IRanges(starts, ends))
  ## Merge columns
  for(coln in colnames(GenomicRanges::mcols)){
    if(coln == 'GT'){
      svs.merged$GT = 'hom'
    } else if(coln == 'QUAL'){ # Should we use the average quality ???
       svs.merged$QUAL = (svs$QUAL[ol$queryHits] + svs$QUAL[ol$subjectHits])/2
    } else if(coln == 'ALT'){
      svs.merged$ALT = svs$ALT[ol$queryHits]
    } else if(coln == 'type'){
      svs.merged$type = svs$type[ol$queryHits]
    } else if(coln == 'size'){
      svs.merged$size = (svs$size[ol$queryHits] + svs$size[ol$subjectHits]) / 2
    } else if(coln == 'ref.cov'){
      svs.merged$ref.cov = svs$ref.cov[ol$queryHits] + svs$ref.cov[ol$subjectHits]
    } else if(coln == 'alt.cov'){
      svs.merged$alt.cov = svs$alt.cov[ol$queryHits] + svs$alt.cov[ol$subjectHits]
    }
  }
  ## Remove merged pairs and add new SVs
  svs = svs[setdiff(1:length(svs), unique(c(ol$queryHits, ol$subjectHits)))]
  svs = c(svs, svs.merged)
  return(svs)
}
