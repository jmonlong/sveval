##' Compute the precision, recall and F1 score using the TP, TP.baseline, FP
##' and FN columns.
##' @title Compute precision, recall and F1 score
##' @param eval.df a data.frame with columns TP, TP.baseline, FP, and FN.
##' @param use.calls.tp.everywhere use the TP count as number of calls matched even for the recall computation
##' @return the input data.frame with 3 new columns precision, recall and F1.
##' @author Jean Monlong
##' @export
prf <- function(eval.df, use.calls.tp.everywhere=TRUE){
  eval.df$precision = eval.df$TP / (eval.df$TP + eval.df$FP)
  eval.df$precision = round(eval.df$precision, 4)
  if(use.calls.tp.everywhere){
    eval.df$recall = eval.df$TP / (eval.df$TP + eval.df$FN)
  } else {
    eval.df$recall = eval.df$TP.baseline / (eval.df$TP.baseline + eval.df$FN)
  }
  eval.df$recall = round(eval.df$recall, 4)
  eval.df$F1 = 2 * eval.df$precision * eval.df$recall /
    (eval.df$precision + eval.df$recall)
  eval.df$F1 = round(eval.df$F1, 4)
  eval.df$F1 = ifelse(eval.df$precision == 0 & eval.df$recall == 0, 0, eval.df$F1)
  return(eval.df)
}
