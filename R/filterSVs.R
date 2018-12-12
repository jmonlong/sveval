##' @title Filter SVs for regions of interest
##' @param sv.gr the input SVs (e.g. read from \code{readSVvcf})
##' @param regions.gr the regions of interest
##' @param ol.prop minimum proportion of sv.gr that must overlap regions.gr.
##' Default is 0.5
##' @return a subset of sv.gr that overlaps regions.gr
##' @author Jean Monlong
##' @export
filterSVs <- function(sv.gr, regions.gr, ol.prop=.5){
  ol.df = GenomicRanges::findOverlaps(sv.gr, regions.gr)
  ol.df =  as.data.frame(ol.df)
  ol.df$qsw = GenomicRanges::width(GenomicRanges::pintersect(sv.gr[ol.df$queryHits], regions.gr[ol.df$subjectHits]))
  ol.df$qw = GenomicRanges::width(sv.gr[ol.df$queryHits])
  sv.gr[unique(ol.df$queryHits[which(ol.df$qsw >= ol.prop*ol.df$qw)])]
}
