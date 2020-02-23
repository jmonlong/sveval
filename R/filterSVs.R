##' @title Filter SVs for size and regions of interest
##' @param sv.gr the input SVs (e.g. read from \code{readSVvcf})
##' @param regions.gr the regions of interest. Ignored if NULL (default).
##' @param ol.prop minimum proportion of sv.gr that must overlap regions.gr.
##' Default is 0.5
##' @param min.size the minimum SV size to be considered. Default 0.
##' @param max.size the maximum SV size to be considered. Default is Inf.
##' @return a subset of sv.gr that overlaps regions.gr or in the specified size range.
##' @author Jean Monlong
##' @export
filterSVs <- function(sv.gr, regions.gr=NULL, ol.prop=.5, min.size=0, max.size=Inf){
  if(min.size>0 | !is.infinite(max.size)){
    sv.gr = sv.gr[which(sv.gr$size>=min.size &
                        sv.gr$size<=max.size)]
  }
  if(!is.null(regions.gr)){
    ol.df = suppressWarnings(GenomicRanges::findOverlaps(sv.gr, regions.gr))
    ol.df =  as.data.frame(ol.df)
    ol.df$qsw = suppressWarnings(GenomicRanges::width(GenomicRanges::pintersect(sv.gr[ol.df$queryHits], regions.gr[ol.df$subjectHits])))
    ol.df$qw = GenomicRanges::width(sv.gr[ol.df$queryHits])
    sv.gr = sv.gr[unique(ol.df$queryHits[which(ol.df$qsw >= ol.prop*ol.df$qw)])]
  }
  return(sv.gr)
  }
