##' @title Annotate SVs from overlap
##' @param ol.l output of an overlap function (olInsertions or olRanges).
##' @param min.qual the minimum QUAL considered for the calls.
##' @param method the method to annotate the overlap. Either 'coverage' (default) for the
##' cumulative coverage (e.g. to deal with fragmented calls); or 'bipartite' for a 1-to-1
##' matching of variants in the calls and truth sets.
##' @return an updated list with a *cov* column added to the calls and truth sets.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
annotateOl <- function(ol.l, min.qual=0, method=c('coverage', 'bipartite')){
  hq.idx = which(ol.l$calls$QUAL >= min.qual)
  if(length(ol.l$truth)>0){
    ol.l$truth$cov = 0
  }
  if(length(ol.l$calls)>0){
    ol.l$calls$cov = 0
  }
  if('ol' %in% names(ol.l)){
    ## filter based on quality
    ol.l$ol = ol.l$ol[which(ol.l$ol$call.idx %in% hq.idx),]
    if(length(ol.l$ol)>0 && nrow(ol.l$ol)>0){
      if(method[1]=='coverage'){
        ## Coverage on truth set
        ins.truth.cov = ol.l$ol %>%
          dplyr::group_by(.data$truth.idx) %>%
          dplyr::summarize(cov=sum(.data$truth.cov))
        ol.l$truth$cov[ins.truth.cov$truth.idx] = ins.truth.cov$cov
        ## Coverage on call set
        ins.call.cov = ol.l$ol %>%
          dplyr::group_by(.data$call.idx) %>%
          dplyr::summarize(cov=sum(.data$call.cov))
        ol.l$calls$cov[ins.call.cov$call.idx] = ins.call.cov$cov
      } else if(method[1] == 'bipartite'){
        ## label both sets
        ol.l$ol$call.lab = paste0('call_', ol.l$ol$call.idx)
        ol.l$ol$truth.lab = paste0('truth_', ol.l$ol$truth.idx)
        ## list edges
        edge.chars = rbind(ol.l$ol$call.lab, ol.l$ol$truth.lab)
        edge.chars = as.character(as.vector(edge.chars))
        ## make graph
        gg = igraph::make_graph(edge.chars, directed=FALSE)
        ## specify group for bipartite clustering
        igraph::vertex_attr(gg)$type = grepl('call_', igraph::vertex_attr(gg)$name)
        ## weight as reciprocal overlap
        calls.size = ol.l$calls$size[ol.l$ol$call.idx]
        truth.size = ol.l$truth$size[ol.l$ol$truth.idx]
        igraph::edge_attr(gg)$weight = ifelse(calls.size>truth.size,
                                              truth.size/calls.size,
                                              calls.size/truth.size)
        ## bipartite clustering
        bpm = igraph::max_bipartite_match(gg)$matching
        bpm = bpm[which(!is.na(bpm))]
        matched.labs = c(names(bpm), as.character(bpm))
        ## translate into "coverage" for the rest of the pipeline
        truth.idx.m = ol.l$ol$truth.idx[which(ol.l$ol$truth.lab %in% matched.labs)]
        ol.l$truth$cov[truth.idx.m] = ol.l$truth$size[truth.idx.m]
        call.idx.m = ol.l$ol$call.idx[which(ol.l$ol$call.lab %in% matched.labs)]
        ol.l$calls$cov[call.idx.m] = ol.l$calls$size[call.idx.m]
      } else {
        stop('"method=" must be one of "coverage" or "bipartite".')
      }
    }
  } else if('rol.gr' %in% names(ol.l)){
    rol.gr = ol.l$rol.gr[which(ol.l$rol.gr$call.idx %in% hq.idx)]
    if(length(rol.gr)>0){
      if(method[1]=='coverage'){
        ## Overlap coverage on the truth set
        gr.l = GenomicRanges::GRangesList(GenomicRanges::split(rol.gr, rol.gr$truth.idx))
        gr.l = GenomicRanges::reduce(gr.l)
        cov.l = lapply(GenomicRanges::width(gr.l), sum)
        ol.l$truth$cov[as.numeric(names(cov.l))] = unlist(cov.l)
        ## Overlap coverge on the call set
        gr.l = GenomicRanges::GRangesList(GenomicRanges::split(rol.gr, rol.gr$call.idx))
        gr.l = GenomicRanges::reduce(gr.l)
        cov.l = lapply(GenomicRanges::width(gr.l), sum)
        ol.l$calls$cov[as.numeric(names(cov.l))] = unlist(cov.l)
      } else if(method[1] == 'bipartite'){
        ## label both sets
        rol.gr$call.lab = paste0('call_', rol.gr$call.idx)
        rol.gr$truth.lab = paste0('truth_', rol.gr$truth.idx)
        ## list edges
        edge.chars = rbind(rol.gr$call.lab, rol.gr$truth.lab)
        edge.chars = as.character(as.vector(edge.chars))
        ## make graph
        gg = igraph::make_graph(edge.chars, directed=FALSE)
        ## specify group for bipartite clustering
        igraph::vertex_attr(gg)$type = grepl('call_', igraph::vertex_attr(gg)$name)
        ## weight as reciprocal overlap
        calls.size = ol.l$calls$size[rol.gr$call.idx]
        truth.size = ol.l$truth$size[rol.gr$truth.idx]
        ol.size = GenomicRanges::width(rol.gr)
        igraph::edge_attr(gg)$weight = ifelse(calls.size>truth.size,
                                              ol.size/calls.size,
                                              ol.size/truth.size)
        ## bipartite clustering
        bpm = igraph::max_bipartite_match(gg)$matching
        bpm = bpm[which(!is.na(bpm))]
        matched.labs = c(names(bpm), as.character(bpm))
        ## translate into "coverage" for the rest of the pipeline
        truth.idx.m = rol.gr$truth.idx[which(rol.gr$truth.lab %in% matched.labs)]
        ol.l$truth$cov[truth.idx.m] = ol.l$truth$size[truth.idx.m]
        call.idx.m = rol.gr$call.idx[which(rol.gr$call.lab %in% matched.labs)]
        ol.l$calls$cov[call.idx.m] = ol.l$calls$size[call.idx.m]
      } else {
        stop('"method=" must be one of "coverage" or "bipartite".')
      }
    }
  }
  ol.l$calls = ol.l$calls[hq.idx]
  return(ol.l)
}
