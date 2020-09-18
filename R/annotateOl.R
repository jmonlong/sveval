##' @title Annotate SVs from overlap
##' @param ol.l output of an overlap function (olInsertions or olRanges).
##' @param min.qual the minimum quality considered for the calls.
##' @param min.cov the minimum coverage to be considered a match. Default is 0.5
##' @param method the method to annotate the overlap. Either 'coverage' (default) for the
##' cumulative coverage (e.g. to deal with fragmented calls); or 'bipartite' for a 1-to-1
##' matching of variants in the calls and truth sets.
##' @return an updated list with a *cov* column added to the calls and truth sets.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
annotateOl <- function(ol.l, min.qual=0, min.cov=.5, method=c('coverage', 'bipartite')){
  hq.idx = which(ol.l$calls$qual >= min.qual)
  if(length(ol.l$truth)>0){
    ol.l$truth$ol = FALSE
  }
  if(length(ol.l$calls)>0){
    ol.l$calls$ol = FALSE
  }
  if('ol' %in% names(ol.l)){
    ## this is for insertions where we have a list of pairs of insertions that cluster together
    ## filter based on quality
    ol.l$ol = ol.l$ol[which(ol.l$ol$call.idx %in% hq.idx),]
    if(length(ol.l$ol)>0 && nrow(ol.l$ol)>0){
      if(method[1]=='coverage'){
        ## Coverage on truth set
        ins.truth.cov = ol.l$ol %>%
          dplyr::group_by(.data$truth.idx) %>%
          dplyr::summarize(cov=sum(.data$truth.cov))
        cov.temp = rep(0, length(ol.l$truth))
        cov.temp[ins.truth.cov$truth.idx] = ins.truth.cov$cov
        ol.l$truth$ol = cov.temp / ol.l$truth$size >= min.cov
        ## Coverage on call set
        ins.call.cov = ol.l$ol %>%
          dplyr::group_by(.data$call.idx) %>%
          dplyr::summarize(cov=sum(.data$call.cov))
        cov.temp = rep(0, length(ol.l$calls))
        cov.temp[ins.call.cov$call.idx] = ins.call.cov$cov
        ol.l$calls$ol = cov.temp / ol.l$calls$size >= min.cov
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
        ## components with only 2 variants are matched, no need to use bipartite
        comp.m = igraph::components(gg)$membership
        comp.mt = table(comp.m)
        comp2 = names(comp.mt)[which(comp.mt==2)]
        comp2.v = names(comp.m)[comp.m %in% comp2]
        gg = igraph::delete_vertices(gg, comp2.v)
        ## bipartite clustering on variants in larger components
        bpm = igraph::max_bipartite_match(gg)$matching
        bpm = bpm[which(!is.na(bpm))]
        ## matched variants: variants in 2-variants comp or bipartite matched
        matched.labs = c(comp2.v, names(bpm), as.character(bpm))
        ## annotate the 'ol' column in each SV set
        truth.idx.m = ol.l$ol$truth.idx[which(ol.l$ol$truth.lab %in% matched.labs)]
        ol.l$truth$ol[truth.idx.m] = TRUE
        call.idx.m = ol.l$ol$call.idx[which(ol.l$ol$call.lab %in% matched.labs)]
        ol.l$calls$ol[call.idx.m] = TRUE
      } else {
        stop('"method=" must be one of "coverage" or "bipartite".')
      }
    } else {
      ol.o$ol = data.frame()
    }
  } else if('rol.gr' %in% names(ol.l)){
    ## this is when comparing ranges, e.g. deletions or inversions
    rol.gr = ol.l$rol.gr[which(ol.l$rol.gr$call.idx %in% hq.idx)]
    if(length(rol.gr)>0){
      if(method[1]=='coverage'){
        ## Overlap coverage on the truth set
        gr.l = GenomicRanges::GRangesList(GenomicRanges::split(rol.gr, rol.gr$truth.idx))
        gr.l = GenomicRanges::reduce(gr.l)
        cov.l = lapply(GenomicRanges::width(gr.l), sum)
        cov.temp = rep(0, length(ol.l$truth))
        cov.temp[as.numeric(names(cov.l))] = unlist(cov.l)
        ol.l$truth$ol = cov.temp /ol.l$truth$size >= min.cov
        ## Overlap coverge on the call set
        gr.l = GenomicRanges::GRangesList(GenomicRanges::split(rol.gr, rol.gr$call.idx))
        gr.l = GenomicRanges::reduce(gr.l)
        cov.l = lapply(GenomicRanges::width(gr.l), sum)
        cov.temp = rep(0, length(ol.l$calls))
        cov.temp[as.numeric(names(cov.l))] = unlist(cov.l)
        ol.l$calls$ol = cov.temp /ol.l$calls$size >= min.cov
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
        ## components with only 2 variants are matched, no need to use bipartite
        comp.m = igraph::components(gg)$membership
        comp.mt = table(comp.m)
        comp2 = names(comp.mt)[which(comp.mt==2)]
        comp2.v = names(comp.m)[comp.m %in% comp2]
        gg = igraph::delete_vertices(gg, comp2.v)
        ## bipartite clustering on variants in larger components
        bpm = igraph::max_bipartite_match(gg)$matching
        bpm = bpm[which(!is.na(bpm))]
        ## matched variants: variants in 2-variants comp or bipartite matched
        matched.labs = c(comp2.v, names(bpm), as.character(bpm))
        ## annotate the 'ol' column in each SV set
        truth.idx.m = rol.gr$truth.idx[which(rol.gr$truth.lab %in% matched.labs)]
        ol.l$truth$ol[truth.idx.m] = TRUE
        call.idx.m = rol.gr$call.idx[which(rol.gr$call.lab %in% matched.labs)]
        ol.l$calls$ol[call.idx.m] = TRUE
      } else {
        stop('"method=" must be one of "coverage" or "bipartite".')
      }
    }
  }
  ol.l$calls = ol.l$calls[hq.idx]
  return(ol.l)
}
