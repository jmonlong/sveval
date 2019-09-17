#' Evaluate SV in a call set against a truth set using overlap-based approaches and sequence comparison for insertions.
#'
#' \tabular{ll}{
#' Package: \tab sveval\cr
#' Type: \tab Package\cr
#' Version: \tab 1.2.2\cr
#' Date: \tab 2019-09-16\cr
#' License: \tab MIT\cr
#' }
#' @docType package
#' @name sveval-package
#' @title SV evaluation
#' @author Jean Monlong \email{jmonlong@ucsc.edu}
#' @seealso \url{http://www.github.com/jmonlong/sveval}
##' @examples
##' \dontrun{
##' eval = svevalOl('calls.vcf', 'truth.vcf')
##' plot_prcurve(eval$curve)
##'
##' # Comparing multiple methods
##' eval.1 = svevalOl('calls1.vcf', 'truth.vcf')
##' eval.2 = svevalOl('calls2.vcf', 'truth.vcf')
##' plot_prcurve(list(eval.1$curve, eval.2$curve), labels=c('method1', 'method2'))
##' }
NULL
