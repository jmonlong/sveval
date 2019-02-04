##' @title Annotate SVs from overlap
##' @param ol.l output of an overlap function (olInsertions or olRanges).
##' @param min.qual the minimum QUAL considered for the calls.
##' @return an updated list with a *cov* column added to the calls and truth sets.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
annotateOl <- function(ol.l, min.qual=0){
  hq.idx = which(ol.l$calls$QUAL >= min.qual)
  if(length(ol.l$truth)>0){
    ol.l$truth$cov = 0
  }
  if(length(ol.l$calls)>0){
    ol.l$calls$cov = 0
  }
  if('ol' %in% names(ol.l)){
    ## Coverage on truth set
    if(length(ol.l$truth)>0 & length(ol.l$calls)>0){
      ins.truth.cov = ol.l$ol %>%
        dplyr::filter(.data$call.idx %in% hq.idx) %>% 
        dplyr::group_by(.data$truth.idx) %>%
        dplyr::summarize(cov=sum(.data$truth.cov))
      ol.l$truth$cov[ins.truth.cov$truth.idx] = ins.truth.cov$cov
      ## Coverage on call set
      ol.l$calls$cov = 0
      ins.call.cov = ol.l$ol %>%
        dplyr::filter(.data$call.idx %in% hq.idx) %>% 
        dplyr::group_by(.data$call.idx) %>%
        dplyr::summarize(cov=sum(.data$call.cov))
      ol.l$calls$cov[ins.call.cov$call.idx] = ins.call.cov$cov
    }
  } else if('rol.gr' %in% names(ol.l)){
    rol.gr = ol.l$rol.gr[which(ol.l$rol.gr$call.idx %in% hq.idx)]
    ## Overlap coverage on the truth set
    if(length(ol.l$truth)>0 & length(ol.l$calls)>0){
      gr.l = GenomicRanges::GRangesList(GenomicRanges::split(rol.gr, rol.gr$truth.idx))
      gr.l = GenomicRanges::reduce(gr.l)
      cov.l = lapply(GenomicRanges::width(gr.l), sum)
      ol.l$truth$cov[as.numeric(names(cov.l))] = unlist(cov.l)
      ## Overlap coverge on the call set
      gr.l = GenomicRanges::GRangesList(GenomicRanges::split(rol.gr, rol.gr$call.idx))
      gr.l = GenomicRanges::reduce(gr.l)
      cov.l = lapply(GenomicRanges::width(gr.l), sum)
      ol.l$calls$cov[as.numeric(names(cov.l))] = unlist(cov.l)
    }
  }
  ol.l$calls = ol.l$calls[hq.idx]
  return(ol.l)
}
