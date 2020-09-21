##' @title Annotate the overlap between SVs using reciprocal overlap and bipartite clustering
##' @param ol.gr overlaps prepared with \code{prepareOl} function.
##' @param min.ol the minimum overlap to be considered a potential match. Default is 0.5
##' @return an updated and filtered version of the input ol.gr
##' @author Jean Monlong
##' @keywords internal
annotateOl.bipartite <- function(ol.gr, min.ol=.5){
  ol.gr$queryOl = ol.gr$subjectOl = FALSE

  ## compute the overlap, here simple reciprocal overlap
  ol.gr$olScore = ifelse(ol.gr$subjectSize > ol.gr$querySize,
                         ol.gr$interSize / ol.gr$subjectSize,
                         ol.gr$interSize / ol.gr$querySize)

  ## remove ovelap shorter than min.ol
  ol.gr = ol.gr[which(ol.gr$olScore >= min.ol)]

  if(length(ol.gr)> 0){
    ## label both sets
    ol.gr$queryLab = paste0('q_', ol.gr$queryHits)
    ol.gr$subjectLab = paste0('s_', ol.gr$subjectHits)
    
    ## list edges
    edge.chars = rbind(ol.gr$queryLab, ol.gr$subjectLab)
    edge.chars = as.character(as.vector(edge.chars))

    ## make graph
    gg = igraph::make_graph(edge.chars, directed=FALSE)
    
    ## specify group for bipartite clustering
    igraph::vertex_attr(gg)$type = grepl('q_', igraph::vertex_attr(gg)$name)

    ## weight as reciprocal overlap
    igraph::edge_attr(gg)$weight = ol.gr$olScore

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
    ol.gr$queryOl[which(ol.gr$queryLab %in% matched.labs)] = TRUE
    ol.gr$subjectOl[which(ol.gr$subjectLab %in% matched.labs)] = TRUE
  }
  
  ## remove temporary columns
  ol.gr$queryLab = ol.gr$subjectLab = NULL
  
  return(ol.gr[which(ol.gr$queryOl & ol.gr$subjectOl)])
}
