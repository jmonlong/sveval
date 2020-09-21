##' Overlaps that were prepared by \code{prepateOl} and then potentially filtered
##' (by size, region of interest, call quality), are filtered to match at least
##' a minimum overlap criteria. Different methods can be used to match variants,
##' see below.
##'
##' The 'method' is either 'coverage' (default) for the
##' cumulative coverage (e.g. to deal with fragmented calls); or 'bipartite' for a 1-to-1
##' matching of variants in the calls and truth sets.
##'
##' For some overlap methods (e.g. reciprocal coverage), a match is not always reciprocal.
##' For example, a smaller deletion can be covered by a larger one while still not
##' covering the large one enough. The means the small deletions is "overlapped" (here covered)
##' but not the large one. To report this kind of overlap, two additional columns are added:
##' 'queryOl' and 'subjectOl'. Hence, during evaluation it is not sufficient to cound SVs in
##' the overlap object but rather the ones with queryOl==TRUE or subjectOl==TRUE. When using
##' other overlap approaches (reciprocal overlap and/or bipartite clustering), both columns should
##' always be TRUE.
##' 
##' @title Annotate the overlap between SVs
##' @param ol.gr overlaps prepared with \code{prepareOl} function.
##' @param min.cov the minimum overlap/coverage to be considered a match. Default is 0.5
##' @param method the method to annotate the overlap. See details.
##' @return an updated and filtered version of the input ol.gr. The potential new columns include:
##' \item{queryOl}{should the query be counted as "overlapped"}
##' \item{subjectOl}{should the subject be counted as "overlapped"}
##' \item{olScore}{the overlap score (usually the value of the reciprocal overlap)}
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
annotateOl <- function(ol.gr, min.ol=.5, method=c('coverage', 'reciprocal', 'bipartite')){
  if(length(ol.gr)==0){
    return(ol.gr)
  }
  if(method[1] == 'coverage'){
    return(annotateOl.coverage(ol.gr, min.cov=min.ol))
  } else if (method[1] == 'bipartite'){
    return(annotateOl.bipartite(ol.gr, min.ol=min.ol))
  } else if (method[1] == 'reciprocal'){
    return(annotateOl.reciprocal(ol.gr, min.ol=min.ol))
  } else {
    stop('Overlap method must be one of "coverage" or "bipartite"')
  }
}
