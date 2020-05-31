##' @title Overlap calls and truth set for insertions
##' @param calls.gr call set
##' @param truth.gr truth set
##' @param max.ins.gap maximum distance for insertions to be clustered.
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @return a list with
##' \item{calls}{the input calls}
##' \item{truth}{the input truth}
##' \item{ol}{a data.frame with the overlap information}
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
olInsertions <- function(calls.gr, truth.gr, max.ins.gap=1, 
                       ins.seq.comp=FALSE, nb.cores=1){
  calls.ins  = calls.gr[which(calls.gr$type=='INS')]
  truth.ins  = truth.gr[which(truth.gr$type=='INS')]
  ol.ins = NULL
  if(length(calls.ins)>0 & length(truth.ins)>0){
    ## Cluster insertions
    ol.ins = GenomicRanges::findOverlaps(truth.ins, calls.ins,
                                         maxgap=max.ins.gap)
    ol.ins = as.data.frame(ol.ins)
    if(ins.seq.comp){
      ## Sequence comparison
      if(!('alt' %in% colnames(GenomicRanges::mcols(truth.ins))) |
         !('alt' %in% colnames(GenomicRanges::mcols(calls.ins)))){
        stop('Missing sequence information. Did you run use "keep.ins.seq" when reading the VCF?')
      }
      truth.seq = truth.ins$alt[ol.ins$queryHits]
      calls.seq = calls.ins$alt[ol.ins$subjectHits]
      if(nb.cores > 1){
        chunk.idx = tapply(1:length(truth.seq),
                           cut(1:length(truth.seq), nb.cores),
                           identity,
                           simplify=FALSE)
        res = parallel::mclapply(chunk.idx, function(ii){
          pas = Biostrings::pairwiseAlignment(truth.seq[ii], calls.seq[ii],
                                              type='local')
          ## Biostrings::nchar(pas)
          Biostrings::nmatch(pas)
        }, mc.cores=nb.cores)
        ol.ins$call.cov = ol.ins$truth.cov = unlist(res)
      } else {
        pas = Biostrings::pairwiseAlignment(truth.seq, calls.seq)
        ## ol.ins$cov = Biostrings::nchar(pas)
        ol.ins$call.cov = ol.ins$truth.cov = Biostrings::nmatch(pas)
      }
    } else {
      ## Size comparison
      ol.ins$call.cov = truth.ins$size[ol.ins$queryHits]
      ol.ins$truth.cov = calls.ins$size[ol.ins$subjectHits]
    }
    ol.ins$truth.idx = ol.ins$queryHits
    ol.ins$call.idx = ol.ins$subjectHits
    ol.ins$queryHits = ol.ins$subjectHits = NULL
  }
  return(list(calls=calls.ins, truth=truth.ins, ol=ol.ins))
}
