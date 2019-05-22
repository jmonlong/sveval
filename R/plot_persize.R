##' @title Recall, precision, F1 per SV size
##' @param eval the output of \code{svevalOl}.
##' @param size.breaks a vector how to break the sizes into classes.
##' @param plot should the function return the plot list. Default is TRUE. If FALSE,
##' returns a data.frame.
##' @return a list of ggplot objects if plot=TRUE (default); a data.frame otherwise.
##' @author Jean Monlong
##' @import ggplot2
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
plot_persize <- function(eval, size.breaks=c(50,100,500,1e3,1e4,Inf), plot=TRUE){
  svs = eval$svs
  size.levels = levels(cut(size.breaks, breaks=size.breaks, include.lowest=TRUE))
  
  ## For each SV type
  eval.df = lapply(names(svs), function(svtype){
    svs = svs[[svtype]]

    ## For each class of variant
    df = lapply(c('TP', 'TP.baseline', 'FP', 'FN'), function(metric){
      svs = svs[[metric]]
      if(length(svs) == 0) return(data.frame(type=svtype, metric=metric, size=size.levels,
                                             n=0, stringsAsFactors=FALSE))
      svs$size.class = cut(svs$size, breaks=size.breaks, include.lowest=TRUE)
      
      ## For each size class
      df = lapply(size.levels, function(sl){
        data.frame(type=svtype, metric=metric, size=sl,
                   n=sum(as.character(svs$size.class) == sl, na.rm=TRUE),
                   stringsAsFactors=FALSE)
      })
      do.call(rbind, df)

    })
    do.call(rbind, df)

  })
  eval.df = do.call(rbind, eval.df)
  
  ## Reformat into one row per size class/type with columns TP, FP, etc
  eval.df = tidyr::spread(eval.df, 'metric', 'n', fill=0)
  
  ## Precision, recall and F1
  eval.df$precision = eval.df$TP / (eval.df$TP + eval.df$FP)
  eval.df$precision = round(eval.df$precision, 4)
  eval.df$recall = eval.df$TP.baseline / (eval.df$TP.baseline + eval.df$FN)
  eval.df$recall = round(eval.df$recall, 4)
  eval.df$F1 = 2 * eval.df$precision * eval.df$recall /
    (eval.df$precision + eval.df$recall)
  eval.df$F1 = round(eval.df$F1, 4)

  if(plot){
    ggp.l = list()
    recall = precision = F1 = type = size = NULL
    ggp.l$f1 = ggplot(eval.df, aes(x=size, y=F1, color=type)) +
      geom_line() + geom_point() + theme_bw()
    ggp.l$precision = ggplot(eval.df, aes(x=size, y=precision, color=type)) +
      geom_line() + geom_point() + theme_bw()
    ggp.l$recall = ggplot(eval.df, aes(x=size, y=recall, color=type)) +
      geom_line() + geom_point() + theme_bw()
    return(ggp.l)
  } else {
    return(eval.df)
  } 
}
