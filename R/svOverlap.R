##' @title Overlap and annotate SV sets with coverage metrics
##' @param query a GRanges object with SVs
##' @param subject another GRanges object with SVs
##' @param max.ins.dist maximum distance for insertions to be clustered. Default is 20.
##' @param min.cov the minimum coverage to be considered a match. Default is 0.5
##' @param min.del.rol minimum reciprocal overlap for deletions. Default is 0.1
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @return a list with:
##' \item{query}{the query GRanges object annotated}
##' \item{subject}{the subject GRanges object annotated}
##' @author Jean Monlong
##' @export
svOverlap <- function(query, subject, max.ins.dist=20, min.cov=.5,
                     min.del.rol=.1, ins.seq.comp=FALSE, nb.cores=1){

  ## Using query as "truth" and subject as "calls"
  ## in the internal functions
  
  ## Overlap insertions
  ol.ins = suppressWarnings(
    olInsertions(subject, query, max.ins.gap=max.ins.dist,
                 ins.seq.comp=ins.seq.comp, nb.cores=nb.cores)
  )
  ## Insertion annotation for each genotype
  ins.a = annotateOl(ol.ins)

  ## Overlap deletions
  ol.del = suppressWarnings(
    olRanges(subject, query, min.rol=min.del.rol, type='DEL')
  )
  ## Deletion annotation for each genotype
  del.a = annotateOl(ol.del)

  ## Overlap inversions
  ol.inv = suppressWarnings(
    olRanges(subject, query, min.rol=min.del.rol, type='INV')
  )
  ## Inversion annotation for each genotype
  inv.a = annotateOl(ol.inv)

  ## Merge annotated calls/truth sets
  ol.l = c(ins.a, del.a, inv.a)
  calls=do.call(c, lapply(ol.l, function(ll) ll$calls))
  truth=do.call(c, lapply(ol.l, function(ll) ll$truth))

  calls$cov.prop = calls$cov / calls$size
  truth$cov.prop = truth$cov / calls$size

  return(list(query=truth, subject=calls))
}

