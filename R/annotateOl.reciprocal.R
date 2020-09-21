##' @title Annotate the overlap between SVs using the reciprocal overlap method
##' @param ol.gr overlaps prepared with \code{prepareOl} function.
##' @param min.ol the minimum overlap to be considered a potential match. Default is 0.5
##' @return an updated and filtered version of the input ol.gr
##' @author Jean Monlong
##' @keywords internal
annotateOl.reciprocal <- function(ol.gr, min.ol=.5){
  ol.gr$queryOl = ol.gr$subjectOl = FALSE

  ## compute the overlap, here simple reciprocal overlap
  ol.gr$olScore = ifelse(ol.gr$subjectSize > ol.gr$querySize,
                         ol.gr$interSize / ol.gr$subjectSize,
                         ol.gr$interSize / ol.gr$querySize)

  ## remove ovelap shorter than min.ol
  ol.gr$queryOl = ol.gr$subjectOl = ol.gr$olScore >= min.ol
  return(ol.gr[which(ol.gr$queryOl)])
}
