##' A simple ggplot2 representation of variants in a region.
##' The beginning of the variant is represented as a point (shape=SV type). The point is
##' annotated with the variant size. A line outlines the range (e.g. for deletions or
##' inversions). 
##' @title Plot variants in a region
##' @param gr.l a list of GRanges. If named, the names are used as colors in the graph.
##' @param region.gr the region of interest. If NULL (default), all variants are displayed.
##' @param pt.size the point (and line) sizes. Default is 2.
##' @param lab.size the label size. Default is 4
##' @param maxgap the maximum gap allowed when filtering variants in regions. Default is 20.
##' @param scale.legend the size of the scale legend at the bottom. 0 to switch off. Default is 'auto'.
##' @param show.svids should the SV ids be shown on the y-axis. Default is TRUE
##' @param gr.l.2 (Optional) a second list of GRanges, like gr.l but from another run for example. If provided the graph will show a second panel for this run (change names with \code{run.names=}.
##' @param run.names names to use for the panels when a second input, \code{gr.l.2}, is provided. Vector of two characters.
##' @return a ggplot2 object
##' @author Jean Monlong
##' @export
plot_ranges <- function(gr.l, region.gr=NULL, pt.size=2, lab.size=4, maxgap=20, scale.legend='auto', show.svids=TRUE, gr.l.2=NULL, run.names=c('run1', 'run2')){

  ## add number as names if missing
  if(is.null(names(gr.l))){
    names(gr.l) = 1:length(gr.l)
  }
  if(!is.null(gr.l.2) && is.null(names(gr.l.2))){
    names(gr.l.2) = 1:length(gr.l.2)
  }

  ## zoom to region of interest
  if(!is.null(region.gr)){
    gr.l = lapply(gr.l, function(gr){
      IRanges::subsetByOverlaps(gr, region.gr, maxgap=maxgap)
    })
    if(!is.null(gr.l.2)){
      gr.l.2 = lapply(gr.l.2, function(gr){
        IRanges::subsetByOverlaps(gr, region.gr, maxgap=maxgap)
      })
    }
  }
  ## keep elements with variants
  gr.l = gr.l[which(unlist(lapply(gr.l, length))>0)]
  if(!is.null(gr.l.2)){
    gr.l.2 = gr.l.2[which(unlist(lapply(gr.l.2, length))>0)]  
  }
  
  ## make a data.frame
  if(length(gr.l) == 0){
    df = NULL
  } else {
    df = lapply(1:length(gr.l), function(ii){
      dff = GenomicRanges::as.data.frame(gr.l[[ii]])
      dff$set=names(gr.l)[ii]
      if(!('ac' %in% colnames(dff))){
        dff$ac = NA
      }
      if(!('svid' %in% colnames(dff))){
        dff$svid = NA
      }
      dff[, c('seqnames', 'start', 'end', 'type', 'size', 'set', 'ac', 'svid')]
    })
    df = do.call(rbind, df)
  }
  
  ## if necessary, do the same for the second set
  if(!is.null(gr.l.2)){
    if(length(gr.l.2) == 0){
      df2 = NULL
    } else {
      df2 = lapply(1:length(gr.l.2), function(ii){
        dff = GenomicRanges::as.data.frame(gr.l.2[[ii]])
        dff$set=names(gr.l.2)[ii]
        if(!('ac' %in% colnames(dff))){
          dff$ac = NA
        }
        if(!('svid' %in% colnames(dff))){
          dff$svid = NA
        }
        dff[, c('seqnames', 'start', 'end', 'type', 'size', 'set', 'ac', 'svid')]
      })
      df2 = do.call(rbind, df2)
      df2$run = run.names[2]
    }
    if(!is.null(df)){
      df$run = run.names[1]
      df = rbind(df, df2)
    }
  }
  
  if(!is.null(df)){
    ## order by start position and add IDs
    if('run' %in% colnames(df)){
      df = df[order(df$run, df$set, df$ac, df$start, df$size),]
    } else {
      df = df[order(df$set, df$ac, df$start, df$size),]
    }
    if(all(is.na(df$svid))){
      df$svid = 1:nrow(df)
    } else {
      df$svid = factor(df$svid, levels=unique(df$svid))
    }

    ## label: "SIZE AC" or just "SIZE"
    if(any(!is.na(df$ac))){
      df$label = paste0(format(df$size, big.mark=','), 'bp ', df$ac) 
    } else{
      df$label = paste0(format(df$size, big.mark=','), 'bp')
    }
  } else {
    warning("No variants to show here.")
    return(NULL)
  }
  
  ## graph
  start = end = type = svid = set = label = size = NULL
  min.pos = min(df$start, na.rm=TRUE)
  max.pos = max(df$end, na.rm=TRUE)
  pad.pos = max((max.pos-min.pos) * .1, 30)
  ggp = ggplot2::ggplot(df, ggplot2::aes(color=set))
  if(pt.size=='auto'){
    ggp = ggp + ggplot2::geom_point(ggplot2::aes(x=start, y=svid, shape=type, size=2*size)) +
      ggplot2::geom_segment(ggplot2::aes(x=start, y=svid, xend=end, yend=svid, linewidth=size))
  } else {
    ggp = ggp + ggplot2::geom_point(ggplot2::aes(x=start, y=svid, shape=type), size=pt.size*2) +
      ggplot2::geom_segment(ggplot2::aes(x=start, y=svid, xend=end, yend=svid), linewidth=pt.size)
  }
  ggp = ggp + ggplot2::geom_label(ggplot2::aes(x=end, y=svid, label=label), size=lab.size,
                                  hjust=0, vjust=1, show.legend=FALSE) + 
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y=ggplot2::element_blank()) + 
    ggplot2::scale_x_continuous(lim=c(min.pos-pad.pos, max.pos+pad.pos)) +
    ## ggplot2::scale_y_continuous(lim=c(0, max(df$svid)+1)) + 
    ggplot2::xlab(paste('position on', df$seqnames[1], '(bp)')) +
    ggplot2::ylab('variant')

  if(!show.svids){
    ggp = ggp + ggplot2::theme(axis.text.y=ggplot2::element_blank())
  }

  if(scale.legend=='auto'){
    scale.legend=50
    for(xx in c(100, 1e3, 1e4, 1e5, 1e6)){
      if(max.pos-min.pos > xx*10){
        scale.legend = xx
      }
    }
  }
  
  if(scale.legend>0){
    mid.pos = mean(c(min.pos, max.pos))
    ggp = ggp + ggplot2::annotate('text', x=mid.pos, y=0, label=paste(scale.legend, 'bp'), vjust=0) +
      ggplot2::annotate('errorbarh', xmin=mid.pos-scale.legend/2, xmax=mid.pos+scale.legend/2, y=0)
  }

  if('run' %in% colnames(df)){
    ggp = ggp + ggplot2::facet_grid(run~., scales='free')
  }
  
  return(ggp)
}
