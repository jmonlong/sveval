##' @title Annotate the overlap between SVs using the reciprocal coverage method
##' @param ol.gr overlaps prepared with \code{prepareOl} function.
##' @param min.cov the minimum coverage to be considered a match. Default is 0.5
##' @return an updated and filtered version of the input ol.gr
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
annotateOl.coverage <- function(ol.gr, min.cov=.5){
  ol.gr$queryOl = ol.gr$subjectOl = FALSE
  qs.ol.l = lapply(unique(ol.gr$type), function(svtype){
    ol.gr = ol.gr[which(ol.gr$type==svtype)]
    if(svtype=='INS'){
      ol.df = as.data.frame(ol.gr)
      ## Coverage on queries
      query.cov = ol.df %>%
        dplyr::group_by(.data$queryHits, .data$querySize) %>%
        dplyr::summarize(ol.size=sum(.data$subjectSize)) %>%
        dplyr::filter(.data$ol.size / .data$querySize >= min.cov)
      ## Coverage on call set
      subject.cov = ol.df %>%
        dplyr::group_by(.data$subjectHits, .data$subjectSize) %>%
        dplyr::summarize(ol.size=sum(.data$querySize)) %>%
        dplyr::filter(.data$ol.size / .data$subjectSize >= min.cov)
      return(list(query.ol=unique(query.cov$queryHits),
                  subject.ol=unique(subject.cov$subjectHits)))
    } else { ## "ranges" SVs
      ## The intersection ranges in ol.gr need to be "reduced" to avoid double-counting bases.
      ## this is done after grouping the ranges by query or subject id
      ## Overlap coverage on the truth set
      gr.l = GenomicRanges::GRangesList(GenomicRanges::split(ol.gr, ol.gr$queryHits))
      gr.l = GenomicRanges::reduce(gr.l)
      cov.l = lapply(GenomicRanges::width(gr.l), sum)
      query.cov.df = data.frame(queryHits=as.numeric(names(cov.l)), cov=unlist(cov.l))
      query.cov.df = merge(query.cov.df, as.data.frame(ol.gr)[, c('queryHits', 'querySize')])
      query.cov = query.cov.df$queryHits[which(query.cov.df$cov/query.cov.df$querySize>=min.cov)]
      ## Overlap coverge on the call set
      gr.l = GenomicRanges::GRangesList(GenomicRanges::split(ol.gr, ol.gr$subjectHits))
      gr.l = GenomicRanges::reduce(gr.l)
      cov.l = lapply(GenomicRanges::width(gr.l), sum)
      subject.cov.df = data.frame(subjectHits=as.numeric(names(cov.l)), cov=unlist(cov.l))
      subject.cov.df = merge(subject.cov.df, as.data.frame(ol.gr)[, c('subjectHits', 'subjectSize')])
      subject.cov = subject.cov.df$subjectHits[which(subject.cov.df$cov/subject.cov.df$subjectSize>=min.cov)]
      return(list(query.ol=unique(query.cov),
                  subject.ol=unique(subject.cov)))
    }
  })

  ## merge the query ids to keep and filter the input ol.gr
  query.ol = unlist(lapply(qs.ol.l, function(ll) ll$query.ol))
  ol.gr$queryOl[which(ol.gr$queryHits %in% query.ol)] = TRUE

  ## same for the subject ids
  subject.ol = unlist(lapply(qs.ol.l, function(ll) ll$subject.ol))
  ol.gr$subjectOl[which(ol.gr$subjectHits %in% subject.ol)] = TRUE

  ## compute the reciprocal overlap between the variants, for info
  ol.gr$olScore = ifelse(ol.gr$subjectSize > ol.gr$querySize,
                         ol.gr$interSize / ol.gr$subjectSize,
                         ol.gr$interSize / ol.gr$querySize)

  return(ol.gr[which(ol.gr$queryOl | ol.gr$subjectOl)])
}
