##' @title Prepare SV overlaps before annotation
##' @param query a query GRanges object
##' @param subject a subject GRanges object
##' @param min.rol minimum reciprocal overlap for deletions and other "ranges" SVs. Default is 0.1
##' @param max.ins.dist maximum distance for insertions to be clustered.
##' @param range.seq.comp compare sequence instead of overlapping deletions/inversion/etc. Default is FALSE.
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
                      range.seq.comp=FALSE,
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
  logging::loginfo(paste('Preparing overlaps. Genotypes: ', gts))

  ## prepare overlas separately for SV types and genotypes
  ## pair SV type x genotype (to help potentially parallelize)
  svtypes.gts = lapply(svtypes, function(svtype){
    res = lapply(gts, function(gt){
      list(svtype=svtype, gt=gt)
    })
    do.call(list, res)
  })
  svtypes.gts = do.call(c, svtypes.gts)
  
  if(length(svtypes.gts) == 0){
    return(NULL)
  }
  ol.l = parallel::mclapply(1:length(svtypes.gts), function(ii){
    svtype = svtypes.gts[[ii]]$svtype
    logging::loginfo(paste('Preparing overlaps. ', svtype))
    gt = svtypes.gts[[ii]]$gt
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
    logging::loginfo(paste('Preparing overlaps. # query: ', length(query.ss),
                           ', # subject: ', length(subject.ss)))
    ## stop if one of the set is empty
    if(length(query.ss)==0 | length(subject.ss)==0){
      return(NULL)
    }
    ## otherwise, overlap/compare
    if(svtype == 'INS'){
      ## for insertions, cluster and compare size/sequence
      ## Cluster insertions
      ol.ins = GenomicRanges::findOverlaps(query.ss, subject.ss,
                                           maxgap=max.ins.dist) %>%
        as.data.frame %>% 
        dplyr::mutate(querySize=query.ss$size[.data$queryHits],
                      subjectSize=subject.ss$size[.data$subjectHits],
                      interSize=ifelse(.data$subjectSize>.data$querySize, .data$querySize, .data$subjectSize)) %>%
        dplyr::filter(.data$interSize >= min.rol * .data$querySize,
                      .data$interSize >= min.rol * .data$subjectSize)
      if(nrow(ol.ins)==0) return(NULL)
      if(ins.seq.comp){
        ## Sequence comparison
        if(!('alt' %in% colnames(GenomicRanges::mcols(query.ss))) |
           !('alt' %in% colnames(GenomicRanges::mcols(subject.ss)))){
          stop('Missing sequence information. Did you run use "keep.ins.seq" when reading the VCF?')
        } else {
          if(any(grepl('...',  query.ss$alt, fixed=TRUE)) |
             any(grepl('...',  subject.ss$alt, fixed=TRUE))){
            stop('Missing sequence information. Did you run use "keep.ins.seq" when reading the VCF?')
          }
        }
        query.seq = query.ss$alt[ol.ins$queryHits]
        subject.seq = subject.ss$alt[ol.ins$subjectHits]

        ## ol.ins$interSize = ifelse(ol.ins$querySize > ol.ins$subjectSize,
        ##                           ol.ins$querySize, ol.ins$subjectSize)
        ol.ins$interSize = ol.ins$querySize - aldist(query.seq, subject.seq)
      } else {
        ## Size comparison, save the smallest size
        ol.ins$interSize = ifelse(ol.ins$querySize < ol.ins$subjectSize,
                                  ol.ins$querySize, ol.ins$subjectSize)
      }

      ## format GRanges
      ol.gr = query.ss[ol.ins$queryHits]
      GenomicRanges::mcols(ol.gr) = ol.ins
    } else if(svtype == 'BND'){
      ## make GRange record for each breakpoint
      query.ss.bks = suppressWarnings(c(
        GenomicRanges::GRanges(query.ss$CHR2, IRanges::IRanges(query.ss$end2, width=1), idx=1:length(query.ss), bk=1),
        GenomicRanges::GRanges(GenomicRanges::seqnames(query.ss),
                               IRanges::IRanges(GenomicRanges::start(query.ss),
                                                GenomicRanges::end(query.ss)),
                               idx=1:length(query.ss), bk=2)))
      subject.ss.bks = suppressWarnings(c(
        GenomicRanges::GRanges(subject.ss$CHR2, IRanges::IRanges(subject.ss$end2, width=1),
                               idx=1:length(subject.ss), bk=1),
        GenomicRanges::GRanges(GenomicRanges::seqnames(subject.ss),
                               IRanges::IRanges(GenomicRanges::start(subject.ss),
                                                GenomicRanges::end(subject.ss)),
                               idx=1:length(subject.ss), bk=2)))
      print(query.ss.bks)
      print(subject.ss.bks)
      print(max.ins.dist)
      ## overlap breakpoints forcing more than one breakpoint pairs to match
      ol.bk = GenomicRanges::findOverlaps(query.ss.bks, subject.ss.bks,
                                          maxgap=max.ins.dist) %>%
        as.data.frame %>% 
        dplyr::mutate(bk=paste(query.ss.bks$bk[.data$queryHits], subject.ss.bks$bk[.data$subjectHits]),
                      queryHits=query.ss.bks$idx[.data$queryHits], subjectHits=subject.ss.bks$idx[.data$subjectHits]) %>% 
        dplyr::group_by(.data$queryHits, .data$subjectHits) %>%
        dplyr::filter(length(unique(.data$bk))>1) %>%
        dplyr::ungroup() %>%
        dplyr::select(.data$queryHits, .data$subjectHits) %>% unique
      if(nrow(ol.bk)==0) return(NULL)
      logging::loginfo(paste('Preparing overlaps. # BND matches: ', nrow(ol.bk)))
      ## format GRanges
      ol.gr = query.ss[ol.bk$queryHits]
      GenomicRanges::mcols(ol.gr) = NULL
      ol.gr$queryHits = ol.bk$queryHits
      ol.gr$subjectHits = ol.bk$subjectHits
      ol.gr$querySize = ol.gr$subjectSize = ol.gr$interSize = 1
    } else {
      ## for "ranges" SVs, reciprocal overlap
      rol.df = GenomicRanges::findOverlaps(query.ss, subject.ss) %>%
        as.data.frame %>%
        dplyr::mutate(querySize=query.ss$size[.data$queryHits],
                      subjectSize=subject.ss$size[.data$subjectHits],
                      interSize=GenomicRanges::width(GenomicRanges::pintersect(query.ss[.data$queryHits],
                                                                               subject.ss[.data$subjectHits])),
                      interSize=ifelse(.data$interSize>.data$querySize, .data$querySize, .data$interSize),
                      interSize=ifelse(.data$interSize>.data$subjectSize, .data$subjectSize, .data$interSize)) %>%
        dplyr::filter(.data$interSize >= min.rol * .data$querySize,
                      .data$interSize >= min.rol * .data$subjectSize)
      if(nrow(rol.df)==0) return(NULL)
      if(range.seq.comp){
        ## Sequence comparison
        if(!('ref' %in% colnames(GenomicRanges::mcols(query.ss))) |
           !('ref' %in% colnames(GenomicRanges::mcols(subject.ss)))){
          stop('Missing sequence information. Did you run use "keep.ref.seq" when reading the VCF?')
        } else {
          if(any(grepl('...',  query.ss$ref, fixed=TRUE)) |
             any(grepl('...',  subject.ss$ref, fixed=TRUE))){
            stop('Missing sequence information. Did you run use "keep.ref.seq" when reading the VCF?')
          }
        }
        query.seq = query.ss$ref[rol.df$queryHits]
        subject.seq = subject.ss$ref[rol.df$subjectHits]
        ## rol.df$interSize = ifelse(rol.df$querySize > rol.df$subjectSize,
        ##                           rol.df$querySize, rol.df$subjectSize)
        rol.df$interSize = rol.df$querySize - aldist(query.seq, subject.seq)
      }
      ol.gr = GenomicRanges::pintersect(query.ss[rol.df$queryHits],
                                        subject.ss[rol.df$subjectHits])
      GenomicRanges::mcols(ol.gr) = rol.df
    }
    
    ## convert the indexes back to those of the whole GRanges
    ol.gr$queryHits = query.ii[ol.gr$queryHits]
    ol.gr$subjectHits = subject.ii[ol.gr$subjectHits]
    ol.gr$type = svtype
    return(ol.gr)
  }, mc.cores=nb.cores)
 
  return(do.call(c, ol.l))
}
