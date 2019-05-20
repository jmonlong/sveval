##' Compare calls with a truth set and identifies which variants from
##' the truth set specifically not called (genotype ./.).
##'
##' Same overlapping strategy as in \code{svevalOl} although here no-calls
##' are kept and there is no splitting by genotype. 
##' @title Find no-calls variants
##' @param calls.gr call set. A GRanges or the path to a VCF file.
##' @param truth.gr truth set. A GRanges or the path to a VCF file.
##' @param max.ins.dist maximum distance for insertions to be clustered. Default is 20.
##' @param min.cov the minimum coverage to be considered a match. Default is 0.5
##' @param min.del.rol minimum reciprocal overlap for deletions. Default is 0.1
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @param sample.name the name of the sample to use if VCF files given as
##' input. If NULL (default), use first sample.
##' @param check.inv should the sequence of MNV be compared to identify inversions. 
##' @return a data.frame with coordinates and variant ids from the truth set
##' corresponding to no-calls.
##' @author Jean Monlong
##' @export
findNocalls <- function(calls.gr, truth.gr, max.ins.dist=20, min.cov=.5,
                        min.del.rol=.1, ins.seq.comp=FALSE, nb.cores=1,
                        sample.name=NULL, check.inv=FALSE){
  ## Read/keep no-calls variants from the call set
  if(is.character(calls.gr) & length(calls.gr)==1){
    calls.gr = readSVvcf(calls.gr, keep.ins.seq=ins.seq.comp,
                         sample.name=sample.name,
                         check.inv=check.inv, nocalls=TRUE)
  }
  calls.gr = calls.gr[which(calls.gr$GT == './.' | calls.gr$GT == '.')]
  if(length(calls.gr) == 0){
    return(data.frame())
  }
  ## Truth set
  if(is.character(truth.gr) & length(truth.gr)==1){
    truth.gr = readSVvcf(truth.gr, keep.ins.seq=ins.seq.comp,
                         sample.name=sample.name,
                         check.inv=check.inv, keep.ids=TRUE)
  }
  if(length(truth.gr) == 0){
    stop("Truth set has no SVs.")
  }
  
  ## Overlap no-calls with truth set
  ol.ins = suppressWarnings(
    olInsertions(calls.gr, truth.gr, max.ins.gap=max.ins.dist,
                 ins.seq.comp=ins.seq.comp, nb.cores=nb.cores)
  )
  ol.del = suppressWarnings(
    olRanges(calls.gr, truth.gr, min.rol=min.del.rol, type='DEL')
  )
  ol.inv = suppressWarnings(
    olRanges(calls.gr, truth.gr, min.rol=min.del.rol, type='INV')
  )

  ## Annotation with the overlap coverage
  ins.a = annotateOl(ol.ins)
  del.a = annotateOl(ol.del)
  inv.a = annotateOl(ol.inv)
  ol.l = list(ins.a, del.a, inv.a)
  ol.l = list(
    calls=do.call(c, lapply(ol.l, function(ll) ll$calls)),
    truth=do.call(c, lapply(ol.l, function(ll) ll$truth))
  )

  ## Extract no-calls variants 
  eval.o = evalOl(ol.l, min.cov=min.cov)
  nocalls.df = lapply(names(eval.o$regions), function(svtype){
    regs = eval.o$regions[[svtype]]
    tps = regs$TP.baseline
    tps = GenomicRanges::resize(tps, fix='center',
                                width=GenomicRanges::width(tps) + 2*max.ins.dist)
    tps$varid = names(tps)
    tps = as.data.frame(tps)
    tps[, c('seqnames','start','end','varid','type')]
  })
  nocalls.df = do.call(rbind, nocalls.df)
  
  return(nocalls.df)
}
