## Overlap variants and annotate them with the proportion of region covered
## by variants in the other set. Adds a column 'cov'.
##' @title Annotate inputs with overlap coverage
##' @param call.gr call set
##' @param truth.gr truth set
##' @param max.ins.gap maximum distance for insertions to be clustered.
##' @param min.del.rol minimum reciprocal overlap for deletions. Default is 0.1
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @return a list with both inputs annotated with the overlap coverage (column 'cov').
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
annotateOl <- function(call.gr, truth.gr, max.ins.gap=1, min.del.rol=0.1,
                       ins.seq.comp=FALSE, nb.cores=1){
  ## Insertions
  call.ins  = call.gr[which(call.gr$type=='INS')]
  truth.ins  = truth.gr[which(truth.gr$type=='INS')]
  if(length(truth.ins)>0){
    truth.ins$cov = 0
  }
  if(length(call.ins)>0){
    call.ins$cov = 0
  }
  if(length(call.ins)>0 & length(truth.ins)>0){
    message('Annotating insertions.')
    ## Cluster insertions
    ol.df = GenomicRanges::findOverlaps(truth.ins, call.ins,
                                        maxgap=max.ins.gap)
    ol.df = as.data.frame(ol.df)
    if(ins.seq.comp){
      ## Sequence comparison
      if(!('ALT' %in% colnames(GenomicRanges::mcols(truth.ins))) |
         !('ALT' %in% colnames(GenomicRanges::mcols(call.ins)))){
        stop('Missing sequence information. Did you run use "keep.ins.seq" when reading the VCF?')
      }
      message('   Pairwise alignment of inserted sequences.')
      truth.seq = lapply(truth.ins$ALT[ol.df$queryHits], '[', 1)
      truth.seq = do.call(c, truth.seq)
      call.seq = lapply(call.ins$ALT[ol.df$subjectHits], '[', 1)
      call.seq = do.call(c, call.seq)
      if(nb.cores > 1){
        chunk.idx = tapply(1:length(truth.seq),
                           cut(1:length(truth.seq), nb.cores),
                           identity,
                           simplify=FALSE)
        res = parallel::mclapply(chunk.idx, function(ii){
          pas = Biostrings::pairwiseAlignment(truth.seq[ii], call.seq[ii],
                                              type='local')
          ## Biostrings::nchar(pas)
          Biostrings::nmatch(pas)
        }, mc.cores=nb.cores)
        ol.df$cov = unlist(res)
      } else {
        pas = Biostrings::pairwiseAlignment(truth.seq, call.seq)
        ## ol.df$cov = Biostrings::nchar(pas)
        ol.df$cov = Biostrings::nmatch(pas)
      }
      message('   Coverage computation.')
      ## Coverage on truth set
      ins.truth.cov = ol.df %>%
        dplyr::group_by(.data$queryHits) %>%
        dplyr::summarize(cov=sum(.data$cov))
      truth.ins$cov[ins.truth.cov$queryHits] = ins.truth.cov$cov
      ## Coverage on call set
      ins.call.cov = ol.df %>%
        dplyr::group_by(.data$subjectHits) %>%
        dplyr::summarize(cov=sum(.data$cov))
      call.ins$cov[ins.call.cov$subjectHits] = ins.call.cov$cov
    } else {
      message('   Size comparison.')
      ## Size comparison
      ol.df$truth.w = truth.ins$size[ol.df$queryHits]
      ol.df$call.w = call.ins$size[ol.df$subjectHits]
      ## Coverage on truth set
      ins.truth.cov = ol.df %>%
        dplyr::group_by(.data$queryHits) %>%
        dplyr::summarize(call.w=sum(.data$call.w))
      truth.ins$cov[ins.truth.cov$queryHits] = ins.truth.cov$call.w
      ## Coverage on call set
      ins.call.cov = ol.df %>%
        dplyr::group_by(.data$subjectHits) %>%
        dplyr::summarize(truth.w=sum(.data$truth.w))
      call.ins$cov[ins.call.cov$subjectHits] = ins.call.cov$truth.w
    }
  }
  ## Deletions
  call.del  = call.gr[which(call.gr$type=='DEL')]
  truth.del  = truth.gr[which(truth.gr$type=='DEL')]
  if(length(call.del)>0){
    call.del$cov = 0
  }
  if(length(truth.del)>0){
    truth.del$cov = 0
  }
  if(length(call.del)>0 & length(truth.del)>0){
    message('Annotating deletions.')
    ## Select overlap with minimum reciprocal overlap
    message('   Reciprocal overlap...')
    rol.df = GenomicRanges::findOverlaps(truth.del, call.del) %>%
      as.data.frame %>%
      dplyr::mutate(q.w=GenomicRanges::width(truth.del)[.data$queryHits],
                    s.w=GenomicRanges::width(call.del)[.data$subjectHits],
                    ol.w=GenomicRanges::width(GenomicRanges::pintersect(truth.del[.data$queryHits], call.del[.data$subjectHits]))) %>%
      dplyr::filter(.data$ol.w >= min.del.rol * .data$q.w,
                    .data$ol.w >= min.del.rol * .data$s.w)
    rol.gr = GenomicRanges::pintersect(truth.del[rol.df$queryHits],
                                       call.del[rol.df$subjectHits])
    message('   Coverage computation...')
    ## Overlap coverage on the truth set
    gr.l = GenomicRanges::GRangesList(GenomicRanges::split(rol.gr, rol.df$queryHits))
    gr.l = GenomicRanges::reduce(gr.l)
    cov.l = lapply(GenomicRanges::width(gr.l), sum)
    truth.del$cov[as.numeric(names(cov.l))] = unlist(cov.l)
    ## Overlap coverge on the call set
    gr.l = GenomicRanges::GRangesList(GenomicRanges::split(rol.gr, rol.df$subjectHits))
    gr.l = GenomicRanges::reduce(gr.l)
    cov.l = lapply(GenomicRanges::width(gr.l), sum)
    call.del$cov[as.numeric(names(cov.l))] = unlist(cov.l)
  }
  return(list(calls=c(call.ins, call.del), truth=c(truth.ins, truth.del)))
}
