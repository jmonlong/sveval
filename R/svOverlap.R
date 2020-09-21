##' Overlap SVs by SV type using one of the overlap approaches (see below). Variants covering a
##' genomic range (e.g. deletions, duplications, inversions) are overlapped while insertions are
##' clustered (using \code{max.ins.dist}) and their size or sequence (if \code{ins.seq.comp})
##' are compared.
##'
##' Available overlap approaches, passed with \code{method=}, include: reciprocal, coverage,
##' bipartite. If you are using this function directly, you might be interested in the
##' 'reciprocal' method (default). When evaluating SVs versus a truthset, \code{svevalOl} uses
##' the 'coverage' method to compare calls (absence/presence) and 'bipartite' to compare
##' genotypes (when run with the recommended settings).
##'
##' The "reciprocal" method corresponds to the simple reciprocal overlap for the variants covering a
##' genomic range (e.g. deletions, duplications, inversions), or the reciprocal size/sequence
##' similarity for insertions.
##'
##' With the "coverage" approach, a variant needs to be covered enough by variants from the
##' other set to be counted "matched" or "overlapped". Here again, the ranges are overlapped for
##' SV spanning a genomic region while for insertions, the size or aligned sequences are summed.
##'
##' With the "bipartite" approach, the variants are first matched using the reciprocal
##' overlap method (see "reciprocal"), and then matched one-to-one using bipartite clustering.
##' This ensures that a variant in one set is only matched to one variant in the other set.
##' Useful when comparing genotypes for example when redundancy should be penalized.
##' 
##' @title Overlap SVs
##' @param query a GRanges object with SVs
##' @param subject another GRanges object with SVs
##' @param min.ol the minimum overlap/coverage to be considered a match. Default is 0.5
##' @param method the method to annotate the overlap. Either 'coverage' (default) for the
##' cumulative coverage (e.g. to deal with fragmented calls); or 'bipartite' for a 1-to-1
##' matching of variants in the calls and truth sets.
##' @param max.ins.dist maximum distance for insertions to be clustered. Default is 20.
##' @param min.del.rol minimum reciprocal overlap for deletions. Default is 0.1
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @return a GRanges with information about pairs of SVs in query and subject that overlap
##' \item{GRange}{intersected ranges (informative for "ranges" SVs)}
##' \item{queryHits}{the id of the input query}
##' \item{subjectHits}{the id of the input subject}
##' \item{querSize}{the size of the input query}
##' \item{subjectSize}{the size of the input subject}
##' \item{interSize}{the size of the intersection (e.g. range, ins size, ins seq alignment)}
##' \item{olScore}{the overlap score (usually the value of the reciprocal overlap)}
##' \item{type}{the SV type of the pair}
##' @author Jean Monlong
##' @export
svOverlap <- function(query, subject, min.ol=.5,
                      method=c('reciprocal', 'coverage', 'bipartite'),
                      max.ins.dist=20, 
                      min.del.rol=.1, ins.seq.comp=FALSE, nb.cores=1){  
  ## Prepare overlap
  ol.gr = prepareOl(query, subject, min.rol=min.del.rol,
                   max.ins.dist=max.ins.dist,
                   ins.seq.comp=ins.seq.comp, nb.cores=nb.cores)

  ## Annotate the overlaps
  ol.gr = annotateOl(ol.gr, min.ol=min.ol, method=method)
  
  return(ol.gr)
}

