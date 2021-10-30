##' @title Compute evaluation results on a subset of the calls
##' @param eval the output of \code{svevalOl}.
##' @param regions.gr GRanges object with regions of interest. NULL (default) means all SVs.
##' @param accepted.filters vector of the values of the FILTER field to keep. If NULL (default), all values are accepted
##' which means all values are kept
##' @param min.region.ol minimum proportion of variant that must overlap
##' regions.gr. Default is 0.5
##' @param nb.cores number of processors to use. Default is 1.
##' @return a list like for \code{svevalOl}
##' @author Jean Monlong
##' @export
subset_eval <- function(eval, regions.gr=NULL, accepted.filters=NULL, min.region.ol=.5, nb.cores=1){
  if(is.null(regions.gr) & is.null(accepted.filters)){
    return(eval)
  }

  ## list of SVs split by type and then TP/FN/FP 
  svs = eval$svs

  ## filter SVs based on regions or FILTER
  svs.filt = lapply(names(svs), function(svtype){
    svs = svs[[svtype]]
    ## For each class of variant
    res = lapply(names(svs), function(metric){
      svs = svs[[metric]]
      if(!is.null(regions.gr) | !is.null(accepted.filters)){
        ## don't filter with FILTER on the truth SVs
        if(metric %in% c('TP.baseline', 'FN')){
          accepted.filters = NULL
        }
        svs = filterSVs(svs, regions.gr=regions.gr, ol.prop=min.region.ol, accepted.filters=accepted.filters)
      }
      return(svs)
    })
    names(res) = names(svs)
    return(res)
  })
  names(svs.filt) = names(svs)

  ## compute the evaluation for each quality threshold
  qual.ths = unique(eval$curve$qual)
  eval.quals = lapply(names(svs.filt), function(svtype){
    svs = svs.filt[[svtype]]
    ## For each class of variant
    df = parallel::mclapply(qual.ths, function(min.qual){
      data.frame(type=svtype,
                 metric=c('TP', 'FP', 'FN'),
                 n=c(sum(svs$TP$qual>=min.qual),
                     sum(svs$FP$qual>=min.qual),
                     length(svs$FN) + sum(svs$TP$qual<min.qual)),
                 qual=min.qual,
                 stringsAsFactors=FALSE)
    }, mc.cores=nb.cores)
    do.call(rbind, df)
  })
  eval.quals = do.call(rbind, eval.quals)
  ## Reformat into one row per size class/type with columns TP, FP, etc
  eval.quals = tidyr::spread(eval.quals, 'metric', 'n', fill=0)
  ## Precision, recall and F1
  eval.quals = prf(eval.quals, use.calls.tp.everywhere=TRUE)

  ## summary evaluation table
  eval.df = eval.quals[which(eval.quals$qual == eval$mqual.bestf1),]
  eval.df = eval.df[which(eval.df$FN + eval.df$FP + eval.df$TP>0),]

  return(list(eval=eval.df,
              curve=eval.quals,
              svs=svs,
              mqual.bestf1=eval$mqual.bestf1))
}
