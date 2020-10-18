##' Fragmented calls often can't be matched as is to another call set.
##' Calls that are fragmented in multiple pieces could be stitched together if they
##' are close by and with the same genotype and SV type.
##' An homozygous call can also be fragmented in two heterozygous calls. To counter this,
##' we can merge pairs of heterozygous SVs (from the same type) that overlap substantially.
##' @title Stitch and merge heterozygous SVs
##' @param svs.gr input GRanges with SVs, e.g. read from \code{readSVvcf}.
##' @param do.stitch should nearby het SVs of the same type be stitched? Default is TRUE.
##' @param do.merge should similar het SVs of the same type "merged" into one homozygous variant? Default is TRUE.
##' @param stitch.dist the maximum distance allowed for two SVs to be stitched.
##' @param min.rol minimum reciprocal overlap to merge two hets into one hom
##' @param max.ins.dist maximum distance for insertions to be clustered when merging hets.
##' @param range.seq.comp compare sequence instead of overlapping deletions/inversion/etc. Default is FALSE.
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @return a GRanges with het SVs from input svs.gr stitched and/or merged.
##' @author Jean Monlong
stitchMergeHets <- function(svs.gr, do.stitch=TRUE, do.merge=TRUE,
                            stitch.dist=20, min.rol=.8, max.ins.dist=20,
                            range.seq.comp=FALSE, ins.seq.comp=FALSE){
  ## iteratively stitch SVs until no more pairs can be stitched
  iterStitch <- function(svs.gr){
    svs.gr = lapply(unique(svs.gr$type), function(type){
      svs.t = svs.gr[which(svs.gr$type==type)]
      nhets = Inf
      logging::loginfo(paste('Iteratively stitch', type))
      while(length((hets.idx = which(svs.t$ac == 1))) < nhets){
        nhets = length(hets.idx)
        hets = stitchSVs(svs.t[hets.idx], stitch.dist=stitch.dist)
        svs.t = c(hets, svs.t[which(svs.t$ac > 1)])
      }
      return(svs.t)
    })
    do.call(c, svs.gr)
  }
  ## iteratively merge hets until no more pairs can be merged
  iterMerge <- function(svs.gr){
    svs.gr = lapply(unique(svs.gr$type), function(type){
      svs.t = svs.gr[which(svs.gr$type==type)]
      nhets = Inf
      logging::loginfo(paste('Iteratively merge hets', type))
      while(length((hets.idx = which(svs.t$ac == 1))) < nhets){
        nhets = length(hets.idx)
        hets = mergeHets(svs.t[hets.idx], min.rol=min.rol,
                         max.ins.dist=max.ins.dist, range.seq.comp=range.seq.comp, ins.seq.comp=ins.seq.comp)
        svs.t = c(hets, svs.t[which(svs.t$ac > 1)])
      }
      return(svs.t)
    })
    do.call(c, svs.gr)
  }
  ## Stitch hets SVs
  if(do.stitch){
    ## Merge hets once first
    if(do.merge){
      logging::loginfo('Pre-stitching merge hets')
      svs.gr = iterMerge(svs.gr)
    }
    logging::loginfo('Stitch hets')
    svs.gr = iterStitch(svs.gr)
  }
  ## Merge hets
  if(do.merge){
    logging::loginfo('Merge hets')
    svs.gr = iterMerge(svs.gr)
  }
  return(svs.gr)
}
