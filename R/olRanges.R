##' @title Overlap SVs as ranges (for deletions and inversions)
##' @param calls.gr the call set
##' @param truth.gr truth set
##' @param min.rol minimum reciprocal overlap for deletions (and inversions). Default is 0.1
##' @param type the type of SVs to overlap. 'DEL' (default) or 'INV'.
##' @return a list with:
##' \item{calls}{the input calls}
##' \item{truth}{the input truth}
##' \item{rol.gr}{a GRanges with intersection of overlapping regions.}
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
olRanges <- function(calls.gr, truth.gr, min.rol=0.1, type=c('DEL', 'INV')){
  calls.r  = calls.gr[which(calls.gr$type==type)]
  truth.r  = truth.gr[which(truth.gr$type==type)]
  rol.gr = NULL
  if(length(calls.r)>0 & length(truth.r)>0){
    rol.df = GenomicRanges::findOverlaps(truth.r, calls.r) %>%
      as.data.frame %>%
      dplyr::mutate(q.w=GenomicRanges::width(truth.r)[.data$queryHits],
                    s.w=GenomicRanges::width(calls.r)[.data$subjectHits],
                    ol.w=GenomicRanges::width(GenomicRanges::pintersect(truth.r[.data$queryHits], calls.r[.data$subjectHits]))) %>%
      dplyr::filter(.data$ol.w >= min.rol * .data$q.w,
                    .data$ol.w >= min.rol * .data$s.w)
    rol.gr = GenomicRanges::pintersect(truth.r[rol.df$queryHits],
                                       calls.r[rol.df$subjectHits])
    rol.gr$truth.idx = rol.df$queryHits
    rol.gr$call.idx = rol.df$subjectHits
  }
  return(list(calls=calls.r, truth=truth.r, rol.gr=rol.gr))
}
