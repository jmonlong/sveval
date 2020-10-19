##' Equivalent SVs are sometimes recorded as quite different variants because placed at
##' different locations of a short tandem repeat. For example, imagine a large 100 bp
##' tandem repeat in the reference genome. An expansion of 50 bp might be represented
##' as a 50 bp insertion at the beginning of the repeat in the callset but at the end
##' of the repeat in the truth set. Because they are distant by 100 bp they might not
##' match. Instead of increasing the distance threshold too much, this function provides
##' a more flexible way of matching variants by first extending them with nearby simple
##' repeat. In this example, because we know of this tandem repeat, both insertions will
##' be extended to span the full annotated reference repeat, hence ensuring that they are
##' matched and compared (e.g. by reciprocal size or sequence alignment distance)
##' short tandem repeat. 
##' @title Extend SV coordinates using simple repeat annotation
##' @param svs.gr a GRanges object with SV information (e.g. read by \code{readSVvcf})
##' @param simprep.gr a GRanges object with simple repeat information.
##' @param max.dist.sr.join maximum distance to join nearby simple repeats. Default is 5 (bp)
##' @param max.sv.dist maximum distance between SV and simple repeat to apply the extension. Default is 5 (bp)
##' @return an updated GRanges object for svs.gr (start/end extended by simprep.gr)
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
extendSVwithSimpRep <- function(svs.gr, simprep.gr, max.dist.sr.join=5, max.sv.dist=5){
  ## join simple repeats
  simprep.gr = GenomicRanges::reduce(simprep.gr, min.gapwidth=max.dist.sr.join)
  
  ## overlap SVs and simple repeats and compute extended coordinates
  ol.df = GenomicRanges::findOverlaps(svs.gr, simprep.gr) %>%
    as.data.frame %>%
    dplyr::mutate(sv.s=GenomicRanges::start(svs.gr)[.data$queryHits],
                  sr.s=GenomicRanges::start(simprep.gr)[.data$subjectHits],
                  sv.e=GenomicRanges::end(svs.gr)[.data$queryHits],
                  sr.e=GenomicRanges::end(simprep.gr)[.data$subjectHits])
  if(nrow(ol.df)>0){
    ol.df = ol.df %>%
    dplyr::group_by(.data$queryHits) %>%
    dplyr::summarize(start=min(c(.data$sv.s, .data$sr.s)),
                     end=max(c(.data$sv.e, .data$sr.e)))
  }
                  
  ## update coordinates
  GenomicRanges::start(svs.gr)[ol.df$queryHits] = ol.df$start
  GenomicRanges::end(svs.gr)[ol.df$queryHits] = ol.df$end
  
  ## return
  return(svs.gr)
}
