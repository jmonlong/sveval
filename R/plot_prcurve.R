##' @title Create precision-recall graphs
##' @param eval a data.frame, a list of data.frames, or a vector with one or several paths to files with curve information.
##' @param labels the labels when multiple results are to be merged. Same length as the
##' list of vector. Ignored is NULL (default).
##' @return list of ggplot graph objects
##' @author Jean Monlong
##' @import ggplot2
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
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
