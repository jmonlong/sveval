##' Create a precision/recall curve using metrics computed by the \code{svevalOl} function.
##' The \code{svevalOl} function returns a list containing a "curve" data.frame with the
##' evaluation metrics for different quality thresholds.
##'
##' If the input is a data.frame (or list of data.frames) it should be the "curve" element
##' of the list returned by the \code{svevalOl} function. If the input is a character (or
##' a vector of characters), they are considered to be file names and the data will be read
##' from these files.
##'
##' If multiple inputs are given, either using a list of data.frames or a vectors with
##' several filenames, one curve per input will be created. This is to be used to quickly
##' compare several methods. The "labels" parameters can be used to specify a label for
##' each input to use for the graphs.
##' @title Create precision-recall graphs
##' @param eval a data.frame, a list of data.frames, or a vector with one or several paths
##' to files with "curve" information.
##' @param labels the labels to use for each input (when multiple inputs are used). Ignored
##' is NULL (default).
##' @return list of ggplot graph objects
##' @author Jean Monlong
##' @import ggplot2
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
##' @examples
##' \dontrun{
##' eval = svevalOl('calls.vcf', 'truth.vcf')
##' plot_prcurve(eval$curve)
##'
##' # Comparing multiple methods
##' eval.1 = svevalOl('calls1.vcf', 'truth.vcf')
##' eval.2 = svevalOl('calls2.vcf', 'truth.vcf')
##' plot_prcurve(list(eval.1$curve, eval.2$curve), labels=c('method1', 'method2'))
##'
##' # Or if the results were previously written in files
##' plot_prcurve(c('methods1-prcurve.tsv', 'methods2-prcurve.tsv'), labels=c('method1', 'method2'))
##' }
plot_prcurve <- function(eval, labels=NULL){
  recall = precision = type = label = qual = NULL
  ggp.l = list()

  ## Check if eval are file names. If so read files.
  if(is.character(eval)){
    if(length(eval)==1){
      eval = utils::read.table(eval, sep='\t', as.is=TRUE, header=TRUE)
    } else {
      if(is.null(labels)){
        labels = make.names(basename(eval))
      }
      if(length(labels) != length(eval)){
        stop('Length of "eval" and "labels" differ.')
      }
      eval = lapply(1:length(eval), function(ii){
        df = utils::read.table(eval[ii], sep='\t', as.is=TRUE, header=TRUE)
        df$label = labels[ii]
        df
      })
    }
  }

  
  ## Simple graph (one method/label)
  if(is.data.frame(eval)){
    eval = eval[which(!is.na(eval$recall) & !is.na(eval$precision)),]
    ggp.l$pr = eval %>% dplyr::arrange(.data$qual) %>% 
      ggplot(aes(x=recall, y=precision, colour=type)) + geom_path() +
      geom_point(aes(size=qual)) + theme_bw()
  } else {
    ## Check that list of data.frames have a label column or add one
    if(!all(sapply(eval, function(df) 'label' %in% colnames(df)))){
      if(is.null(labels)){
        labels = 1:length(eval)
      }
      eval = lapply(1:length(eval), function(ii){
        df = eval[[ii]]
        df$label = labels[ii]
        df
      })
    }
  
    ## Merging multiple results
    eval = do.call(rbind, eval)
    eval = eval[which(!is.na(eval$recall) & !is.na(eval$precision)),]
    ggp.l = lapply(unique(eval$type), function(svtype){
      eval %>%
        dplyr::filter(.data$type==svtype) %>% 
        dplyr::arrange(.data$qual) %>% 
        ggplot(aes(x=recall, y=precision, colour=label)) + geom_path() +
        geom_point(aes(size=qual)) + theme_bw()
    })
    names(ggp.l) = paste0('pr.', unique(eval$type))    
  }

  return(ggp.l)
}
