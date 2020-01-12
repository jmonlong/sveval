##' A simple ggplot2 representation of variants in a region.
##' The beginning of the variant is represented as a point (shape=SV type). The point is
##' annotated with the variant size. A line outlines the range (e.g. for deletions or
##' inversions). 
##' @title Plot variants in a region
##' @param gr.l a list of GRanges. If named, the names are used to name the graph's panel.
##' @param region.gr the region of interest. If NULL (default), all variants are displayed.
##' @param pt.size the point (and line) sizes. Default is 2.
##' @return a ggplot2 object
##' @author Jean Monlong
##' @export
plot_ranges <- function(gr.l, region.gr=NULL, pt.size=2){
  if(is.null(names(gr.l))){
    names(gr.l) = 1:length(gr.l)
  }
  df = lapply(1:length(gr.l), function(ii){
    dff = GenomicRanges::as.data.frame(gr.l[[ii]])
    dff$set=names(gr.l)[ii]
    dff[, c('seqnames', 'start', 'end', 'type', 'size', 'set')]
  })
  df = do.call(rbind, df)
  df = df[order(df$start),]
  df$id = 1:nrow(df)
  start = end = type = id = set = size = NULL
  ggplot2::ggplot(df, ggplot2::aes(color=set)) +
    ggplot2::geom_point(ggplot2::aes(x=start, y=id, shape=type), size=pt.size*2) +
    ggplot2::geom_segment(ggplot2::aes(x=start, y=id, xend=end, yend=id), size=pt.size) +
    ggplot2::geom_label(ggplot2::aes(x=end, y=id, label=size), hjust=0, vjust=1) + 
    ggplot2::theme_bw() + ggplot2::theme(axis.text.y=ggplot2::element_blank()) + 
    ggplot2::xlab('position (bp)') + ggplot2::ylab('variant')
}
