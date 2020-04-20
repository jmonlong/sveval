##' @title Compute FP, FN, TP
##' @param ol.l list returned by annotation functions like \code{annotateOl}
##' @param min.cov the minimum coverage to be considered a match. Default is 0.5
##' @return a data.frame with FP, FN and TP for each SV type.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
evalOl <- function(ol.l, min.cov=.8){
  if(length(ol.l$calls)==0 | length(ol.l$truth)==0){
    eval.df = data.frame(type=c('Total', 'DEL', 'INS', 'INV'), stringsAsFactors=FALSE)
    eval.df$F1 = eval.df$recall = eval.df$precision = eval.df$FN = eval.df$FP = eval.df$TP.baseline = eval.df$TP = NA
    regs = list()
  } else {
    ## Compute TP and FN from the truth set
    eval.df = GenomicRanges::mcols(ol.l$truth)[,c('type','cov','size')] %>%
      as.data.frame %>% 
      dplyr::group_by(.data$type) %>%
      dplyr::summarize(FN=sum(.data$cov / .data$size < min.cov),
                       TP.baseline=sum(.data$cov / .data$size >= min.cov))
    ## Compute TP and FP from the call set and merge adding NAs if necessary
    eval.df = GenomicRanges::mcols(ol.l$calls)[,c('type','cov','size')] %>%
      as.data.frame %>% 
      dplyr::group_by(.data$type) %>%
      dplyr::summarize(FP=sum(.data$cov / .data$size < min.cov),
                       TP=sum(.data$cov / .data$size >= min.cov)) %>%
      merge(eval.df, all=TRUE)
    ## Convert NAs to 0
    eval.df = as.data.frame(lapply(eval.df, function(x) ifelse(is.na(x), 0, x)))
    ## Keep only some columns
    eval.df = eval.df[, c('type', 'TP', 'TP.baseline', 'FP', 'FN')]
    ## Both types
    eval.df = eval.df %>% dplyr::mutate(type='Total') %>%
      dplyr::group_by(.data$type) %>%
      dplyr::summarize_all(sum) %>% rbind(eval.df)
    ## Precision, recall and F1
    eval.df = prf(eval.df)
    ## Reorder rows to version in toil-vg
    eval.df$type = factor(eval.df$type, levels=c('Total', 'INS', 'DEL', 'INV'))
    eval.df = eval.df[order(eval.df$type),]
    ## FP, TP and FN
    types = unique(ol.l$truth$type)
    regs = lapply(types, function(svtype){
      regs = list()
      ## Truth
      gr = ol.l$truth[which(ol.l$truth$type == svtype)]
      grr = gr[which(gr$cov / gr$size < min.cov)]
      grr$alt = NULL
      regs$FN = grr
      grr = gr[which(gr$cov / gr$size >= min.cov)]
      grr$alt = NULL
      regs$TP.baseline = grr
      ## Calls
      gr = ol.l$calls[which(ol.l$calls$type == svtype)]
      grr = gr[which(gr$cov / gr$size < min.cov)]
      grr$alt = NULL
      regs$FP = grr
      gr$alt = NULL
      grr = gr[which(gr$cov / gr$size >= min.cov)]
      regs$TP = grr
      regs
    })
    names(regs) = types
  }
  return(list(eval=eval.df, regions=regs))
}

