##' Compute the precision, recall and F1 score using the TP, TP.baseline, FP
##' and FN columns.
##' @title Compute precision, recall and F1 score
##' @param eval.df a data.frame with columns TP, TP.baseline, FP, and FN.
##' @return the input data.frame with 3 new columns precision, recall and F1.
##' @author Jean Monlong
##' @export
prf <- function(eval.df){
  eval.df$precision = eval.df$TP / (eval.df$TP + eval.df$FP)
  eval.df$precision = round(eval.df$precision, 4)
  eval.df$recall = eval.df$TP.baseline / (eval.df$TP.baseline + eval.df$FN)
  eval.df$recall = round(eval.df$recall, 4)
  eval.df$F1 = 2 * eval.df$precision * eval.df$recall /
    (eval.df$precision + eval.df$recall)
  eval.df$F1 = round(eval.df$F1, 4)
  eval.df$F1 = ifelse(eval.df$precision == 0 & eval.df$recall == 0, 0, eval.df$F1)
  return(eval.df)
}
