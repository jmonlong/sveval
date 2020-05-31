##' A simple ggplot2 representation of variants in a region.
##' The beginning of the variant is represented as a point (shape=SV type). The point is
##' annotated with the variant size. A line outlines the range (e.g. for deletions or
##' inversions). 
##' @title Plot variants in a region
##' @param gr.l a list of GRanges. If named, the names are used to name the graph's panel.
##' @param region.gr the region of interest. If NULL (default), all variants are displayed.
##' @param pt.size the point (and line) sizes. Default is 2.
##' @param lab.size the label size. Default is 4
##' @param maxgap the maximum gap allowed when filtering variants in regions. Default is 20.
##' @param scale.legend the size of the scale legend at the bottom. 0 to switch off. Default is 50 bp
##' @return a ggplot2 object
##' @author Jean Monlong
##' @export
plot_ranges <- function(gr.l, region.gr=NULL, pt.size=2, lab.size=4, maxgap=20, scale.legend=50){

  ## add number as names if missing
  if(is.null(names(gr.l))){
    names(gr.l) = 1:length(gr.l)
  }

  ## zoom to region of interest
  if(!is.null(region.gr)){
    gr.l = lapply(gr.l, function(gr){
      IRanges::subsetByOverlaps(gr, region.gr, maxgap=maxgap)
    })
  }
  ## keep elements with variants
  gr.l = gr.l[which(unlist(lapply(gr.l, length))>0)]

  ## make a data.frame
  df = lapply(1:length(gr.l), function(ii){
    dff = GenomicRanges::as.data.frame(gr.l[[ii]])
    dff$set=names(gr.l)[ii]
    if(!('ac' %in% colnames(dff))){
      dff$ac = NA
    }
    dff[, c('seqnames', 'start', 'end', 'type', 'size', 'set', 'ac')]
  })
  df = do.call(rbind, df)

  ## order by start position and add IDs
  df = df[order(df$start),]
  df$id = 1:nrow(df)

  ## label: "SIZE AC" or just "SIZE"
  if(any(!is.na(df$ac))){
    df$label = paste(df$size, df$ac) 
  } else{
    df$label = df$size
  }

  ## graph
  start = end = type = id = set = label = NULL
  min.pos = min(df$start, na.rm=TRUE)
  max.pos = max(df$end, na.rm=TRUE)
  pad.pos = max((max.pos-min.pos) * .1, 30)
  ggp = ggplot2::ggplot(df, ggplot2::aes(color=set)) +
    ggplot2::geom_point(ggplot2::aes(x=start, y=id, shape=type), size=pt.size*2) +
    ggplot2::geom_segment(ggplot2::aes(x=start, y=id, xend=end, yend=id), size=pt.size) +
    ggplot2::geom_label(ggplot2::aes(x=end, y=id, label=label), size=lab.size,
                        hjust=0, vjust=1, show.legend=FALSE) + 
    ggplot2::theme_bw() + ggplot2::theme(axis.text.y=ggplot2::element_blank()) +
    ggplot2::scale_x_continuous(lim=c(min.pos-pad.pos, max.pos+pad.pos)) + 
    ggplot2::scale_y_continuous(lim=c(0, max(df$id)+1)) + 
    ggplot2::xlab('position (bp)') + ggplot2::ylab('variant')

  if(scale.legend>0){
    mid.pos = mean(c(min.pos, max.pos))
    ggp = ggp + ggplot2::annotate('text', x=mid.pos, y=0, label=paste(scale.legend, 'bp'), vjust=0) +
      ggplot2::annotate('errorbarh', xmin=mid.pos-scale.legend/2, xmax=mid.pos+scale.legend/2, y=0)
  }
  return(ggp)
}
