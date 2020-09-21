##' @title Filter SVs for size and regions of interest
##' @param sv.gr the input SVs (e.g. read from \code{readSVvcf})
##' @param regions.gr the regions of interest. Ignored if NULL (default).
##' @param ol.prop minimum proportion of sv.gr that must overlap regions.gr.
##' Default is 0.5
##' @param min.size the minimum SV size to be considered. Default 0.
##' @param max.size the maximum SV size to be considered. Default is Inf.
##' @param mark.pass don't actually filter variants but mark the variants with a 'pass' columns. Default is FALSE.
##' @return 
##' \item{if mark.pass=FALSE (default)}{a subset of sv.gr that overlaps regions.gr or in the specified size range.}
##' \item{if mark.pass=TRUE}{the input sv.gr with a new column 'pass' to mark the variants to keep.}
##' @author Jean Monlong
##' @export
filterSVs <- function(sv.gr, regions.gr=NULL, ol.prop=.5, min.size=0, max.size=Inf, mark.pass=FALSE){
  if(length(sv.gr)==0) return(sv.gr)
  ## size selection
  if(min.size>0 | !is.infinite(max.size)){
    pass.idx = which(sv.gr$size>=min.size & sv.gr$size<=max.size)
  } else {
    pass.idx = 1:length(sv.gr)
  }
  ## regions of interest
  if(!is.null(regions.gr)){
    ol.df = suppressWarnings(GenomicRanges::findOverlaps(sv.gr[pass.idx], regions.gr))
    ol.df =  as.data.frame(ol.df)
    ol.df$qsw = suppressWarnings(GenomicRanges::width(GenomicRanges::pintersect(sv.gr[ol.df$queryHits],
                                                                                regions.gr[ol.df$subjectHits])))
    ol.df$qw = GenomicRanges::width(sv.gr[ol.df$queryHits])
    ## sv.gr = sv.gr[unique(ol.df$queryHits[which(ol.df$qsw >= ol.prop*ol.df$qw)])]
    pass.idx = pass.idx[unique(ol.df$queryHits[which(ol.df$qsw / ol.df$qw >= ol.prop)])]
  }
  ## filter or mark variants that passed
  if(mark.pass){
    sv.gr$pass = FALSE
    sv.gr$pass[pass.idx] = TRUE
  } else {
    sv.gr[pass.idx]
  }
  return(sv.gr)
}
