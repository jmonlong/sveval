##' @title Recall, precision, F1 per region
##' @param eval the output of \code{svevalOl}.
##' @param regions.gr GRanges object with regions of interest
##' @param min.region.ol minimum proportion of variant that must overlap
##' regions.gr. Default is 0.5
##' @param plot should the function return the plot list. Default is TRUE. If FALSE,
##' returns a data.frame.
##' @return a list of ggplot objects if plot=TRUE (default); a data.frame otherwise.
##' @author Jean Monlong
##' @import ggplot2
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
plot_perregion <- function(eval, regions.gr, min.region.ol=.5, plot=TRUE){
  svs = eval$svs
  
  ## For each SV type
  eval.df = lapply(names(svs), function(svtype){
    svs = svs[[svtype]]
    ## For each class of variant
    df = lapply(c('TP', 'TP.baseline', 'FP', 'FN'), function(metric){
      svs = svs[[metric]]
      ## Keep variants overlapping regions of interest
      svs = filterSVs(svs, regions.gr=regions.gr, ol.prop=min.region.ol)
      
      data.frame(type=svtype, metric=metric, n=length(svs),
                 stringsAsFactors=FALSE)
    })
    do.call(rbind, df)
  })
  eval.df = do.call(rbind, eval.df)
  
  ## Reformat into one row per size class/type with columns TP, FP, etc
  eval.df = tidyr::spread(eval.df, 'metric', 'n', fill=0)
  
  ## Precision, recall and F1
  eval.df = prf(eval.df)

  if(plot){
    eval.df = eval.df[which(!is.na(eval.df$F1)), ]
    ggp.l = list()
    recall = precision = F1 = type = NULL
    ggp.l$f1 = ggplot(eval.df, aes(x=type, y=F1)) +
      geom_bar(stat='identity') + theme_bw()
    ggp.l$precision = ggplot(eval.df, aes(x=type, y=precision)) +
      geom_bar(stat='identity') + theme_bw()
    ggp.l$recall = ggplot(eval.df, aes(x=type, y=recall)) +
      geom_bar(stat='identity') + theme_bw()
    return(ggp.l)
  } else {
    return(eval.df)
  } 
}
