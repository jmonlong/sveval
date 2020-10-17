##' @title Prepare SV overlaps before annotation
##' @param query a query GRanges object
##' @param subject a subject GRanges object
##' @param min.rol minimum reciprocal overlap for deletions and other "ranges" SVs. Default is 0.1
##' @param max.ins.dist maximum distance for insertions to be clustered.
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @param by.gt should the variants be split by genotype? Default is FALSE, i.e. all variants with
##' an alternate allele (ac>0) is considered 'called'.
##' @return a GRanges with information about pairs of SVs in query and subject that overlap
##' \item{GRange}{intersected ranges (for "ranges" SVs)}
##' \item{queryHits}{the id of the input query}
##' \item{subjectHits}{the id of the input subject}
##' \item{querSize}{the size of the input query}
##' \item{subjectSize}{the size of the input subject}
##' \item{interSize}{the size of the intersection (e.g. range, ins size, ins seq alignment)}
##' \item{type}{the SV type of the pair}
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
prepareOl <- function(query, subject, min.rol=.1, max.ins.dist=1, 
                      ins.seq.comp=FALSE, nb.cores=1, by.gt=FALSE){

  svtypes = intersect(unique(query$type), unique(subject$type))

  ## if no allele count info, assume all SVs should be considered together
  if(is.null(query$ac) & length(query)>0){
    query$ac = 1
  }
  if(!is.null(query$ac) && all(query$ac == -1)){
    query$ac = 1
  }
  if(is.null(subject$ac) & length(subject)>0){
    subject$ac = 1
  }
  if(!is.null(subject$ac) && all(subject$ac == -1)){
    subject$ac = 1
  }

  ## list genotypes to consider
  if(by.gt){
    gts = intersect(query$ac, subject$ac)
  } else {
    gts = 'called'
  }

  ## prepare overlas separately for SV types and genotypes
  ol.l = lapply(svtypes, function(svtype){
    ol.gt.l = lapply(gts, function(gt){
      ## save the indexes of the variants considered now (svtype + GT)
      if(gt == 'called'){
        query.ii = which(query$ac > 0 & query$type == svtype)
        subject.ii = which(subject$ac > 0 & subject$type == svtype)
      } else {
        query.ii = which(query$ac == gt & query$type == svtype)
        subject.ii = which(subject$ac == gt & subject$type == svtype)
      }
      ## subset
      query.ss = query[query.ii]
      subject.ss = subject[subject.ii]
      ## stop if one of the set is empty
      if(length(query.ss)==0 | length(subject.ss)==0){
        return(NULL)
      }
      ## otherwise, overlap
      if(svtype != 'INS'){
        ## for "ranges" SVs, reciprocal overlap
        rol.df = GenomicRanges::findOverlaps(query.ss, subject.ss) %>%
          as.data.frame %>%
          dplyr::mutate(querySize=GenomicRanges::width(query.ss)[.data$queryHits],
                        subjectSize=GenomicRanges::width(subject.ss)[.data$subjectHits],
                        interSize=GenomicRanges::width(GenomicRanges::pintersect(query.ss[.data$queryHits],
                                                                                 subject.ss[.data$subjectHits]))) %>%
          dplyr::filter(.data$interSize >= min.rol * .data$querySize,
                        .data$interSize >= min.rol * .data$subjectSize)
        if(nrow(rol.df)==0) return(NULL)
        ol.gr = GenomicRanges::pintersect(query.ss[rol.df$queryHits],
                                          subject.ss[rol.df$subjectHits])
        GenomicRanges::mcols(ol.gr) = rol.df
      } else {
        ## for insertions, cluster and compare size/sequence
        ## Cluster insertions
        ol.ins = GenomicRanges::findOverlaps(query.ss, subject.ss,
                                             maxgap=max.ins.dist)
        ol.ins = as.data.frame(ol.ins)
        if(nrow(ol.ins)==0) return(NULL)
        ## save the sizes of the insertions
        ol.ins$querySize = query.ss$size[ol.ins$queryHits]
        ol.ins$subjectSize = subject.ss$size[ol.ins$subjectHits]
        if(ins.seq.comp){
          ## Sequence comparison
          if(!('alt' %in% colnames(GenomicRanges::mcols(query.ss))) |
             !('alt' %in% colnames(GenomicRanges::mcols(subject.ss)))){
            stop('Missing sequence information. Did you run use "keep.ins.seq" when reading the VCF?')
          }
          query.seq = query.ss$alt[ol.ins$queryHits]
          subject.seq = subject.ss$alt[ol.ins$subjectHits]
          if(nb.cores > 1){
            chunk.idx = tapply(1:length(query.seq),
                               cut(1:length(query.seq), nb.cores),
                               identity,
                               simplify=FALSE)
            res = parallel::mclapply(chunk.idx, function(ii){
              pas = Biostrings::pairwiseAlignment(query.seq[ii], subject.seq[ii],
                                                  type='local')
              ## Biostrings::nchar(pas)
              Biostrings::nmatch(pas)
            }, mc.cores=nb.cores)
            ol.ins$interSize = unlist(res)
          } else {
            pas = Biostrings::pairwiseAlignment(query.seq, subject.seq)
            ## ol.ins$interSize = Biostrings::nchar(pas)
            ol.ins$interSize = Biostrings::nmatch(pas)
          }
        } else {
          ## Size comparison, save the smallest size
          ol.ins$interSize = ifelse(ol.ins$querySize < ol.ins$subjectSize,
                                    ol.ins$querySize, ol.ins$subjectSize)
        }

        ## format GRanges
        ol.gr = query.ss[ol.ins$queryHits]
        GenomicRanges::mcols(ol.gr) = ol.ins
      }
      
      ## convert the indexes back to those of the whole GRanges
      ol.gr$queryHits = query.ii[ol.gr$queryHits]
      ol.gr$subjectHits = subject.ii[ol.gr$subjectHits]
      return(ol.gr)
    })
    ol.gt.gr = do.call(c, ol.gt.l)
    if(length(ol.gt.gr)==0) return(NULL)
    ol.gt.gr$type = svtype
    return(ol.gt.gr)
  })

  return(do.call(c, ol.l))
}
