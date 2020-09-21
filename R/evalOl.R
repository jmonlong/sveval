##' @title Compute FP, FN, TP
##' @param ol.gr an GRanges with overlap information and filtered/annotated by \code{annotateOl}
##' @param truth.gr the original SVs from the truthset. Correspond to the queryHits in ol.gr. Use
##' only variants with 'pass' column.
##' @param calls.gr the original SVs from the callset. Correspond to the subjectHits in ol.gr. Use
##' only variants with 'pass' column.
##' @return a list with
##' \item{eval}{a data.frame with FP, FN and TP for each SV type.}
##' \item{regions}{a list of FN/FP/TP regions for each SV type}
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
evalOl <- function(ol.gr, truth.gr, calls.gr){
  svtypes = unique(c(truth.gr$type, calls.gr$type))

  ## Add pass column is absent
  if(length(truth.gr)>0 & is.null(truth.gr$pass)) truth.gr$pass = TRUE
  if(length(calls.gr)>0 & is.null(calls.gr$pass)) calls.gr$pass = TRUE
  
  ## Compute TP and FP from the call set
  if(sum(calls.gr$pass)==0){
    eval.call = data.frame(type=c('Total', svtypes), stringsAsFactors=FALSE)
    eval.call$FP = eval.call$TP = 0
  } else {
    calls.gr$ol = FALSE
    if(length(ol.gr)>0){
      calls.gr$ol[ol.gr$subjectHits[which(ol.gr$subjectOl)]] = TRUE
    }
    eval.call = GenomicRanges::mcols(calls.gr)[,c('type', 'pass', 'ol')] %>%
      as.data.frame %>%
      dplyr::filter(.data$pass) %>% 
      dplyr::group_by(.data$type) %>%
      dplyr::summarize(FP=sum(!.data$ol), TP=sum(.data$ol))
    eval.call = eval.call %>% dplyr::mutate(type='Total') %>%
      dplyr::group_by(.data$type) %>%
      dplyr::summarize_all(sum) %>% rbind(eval.call)
  }

  ## Compute TP.baseline and FN from the call set
  if(sum(truth.gr$pass)==0){
    eval.truth = data.frame(type=c('Total', svtypes), stringsAsFactors=FALSE)
    eval.truth$FN = eval.truth$TP.baseline = 0
  } else {
    truth.gr$ol = FALSE
    if(length(ol.gr)>0){
      truth.gr$ol[ol.gr$queryHits[which(ol.gr$queryOl)]] = TRUE
    }
    eval.truth = GenomicRanges::mcols(truth.gr)[,c('type', 'pass', 'ol')] %>%
      as.data.frame %>%
      dplyr::filter(.data$pass) %>% 
      dplyr::group_by(.data$type) %>%
      dplyr::summarize(FN=sum(!.data$ol), TP.baseline=sum(.data$ol))
    eval.truth = eval.truth %>% dplyr::mutate(type='Total') %>%
      dplyr::group_by(.data$type) %>%
      dplyr::summarize_all(sum) %>% rbind(eval.truth)
  }

  ## merge
  eval.df = merge(eval.call, eval.truth, by='type', all=TRUE)
  ## Convert NAs to 0
  eval.df = as.data.frame(lapply(eval.df, function(x) ifelse(is.na(x), 0, x)))
  ## Precision, recall and F1
  eval.df = prf(eval.df)
  ## Reorder rows to version in toil-vg
  eval.df$type = factor(eval.df$type, levels=unique(c('Total', 'INS', 'DEL', 'INV', sort(svtypes))))
  eval.df = eval.df[order(eval.df$type),]
  ## FP, TP and FN
  regs = lapply(svtypes, function(svtype){
    regs = list()
    ## Truth
    gr = truth.gr[which(truth.gr$type == svtype & truth.gr$pass)]
    if(length(gr)==0){
      regs$FN = regs$TP.baseline = gr
    } else {
      grr = gr[which(!gr$ol)]
      grr$alt = NULL
      regs$FN = grr
      grr = gr[which(gr$ol)]
      grr$alt = NULL
      regs$TP.baseline = grr
    }
    ## Calls
    gr = calls.gr[which(calls.gr$type == svtype & calls.gr$pass)]
    if(length(gr)==0){
      regs$FP = regs$TP = gr
    } else {
      grr = gr[which(!gr$ol)]
      grr$alt = NULL
      regs$FP = grr
      gr$alt = NULL
      grr = gr[which(gr$ol)]
      regs$TP = grr
    }
    regs
  })
  names(regs) = svtypes

  return(list(eval=eval.df, regions=regs))
}

