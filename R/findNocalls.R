##' Compare calls with a truth set and identifies which variants from
##' the truth set specifically not called (genotype ./.).
##'
##' Same overlapping strategy as in \code{svevalOl} although here no-calls
##' are kept and there is no splitting by genotype. 
##' @title Find no-calls variants
##' @param calls.gr call set. A GRanges or the path to a VCF file.
##' @param truth.gr truth set. A GRanges or the path to a VCF file.
##' @param max.ins.dist maximum distance for insertions to be clustered. Default is 20.
##' @param min.ol the minimum overlap/coverage to be considered a match. Default is 0.5
##' @param min.del.rol minimum reciprocal overlap for deletions. Default is 0.1
##' @param range.seq.comp compare sequence instead of only overlapping deletions/inversions/etc. Default is FALSE.
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @param sample.name the name of the sample to use if VCF files given as
##' input. If NULL (default), use first sample.
##' @param check.inv should the sequence of MNV be compared to identify inversions. 
##' @return a data.frame with coordinates and variant ids from the truth set
##' corresponding to no-calls.
##' @param method the method to annotate the overlap. Either 'coverage' (default) for the
##' cumulative coverage (e.g. to deal with fragmented calls); or 'bipartite' for a 1-to-1
##' matching of variants in the calls and truth sets.
##' @author Jean Monlong
##' @export
findNocalls <- function(calls.gr, truth.gr, max.ins.dist=20, min.ol=.5,
                        min.del.rol=.1, range.seq.comp=FALSE, ins.seq.comp=FALSE, nb.cores=1,
                        sample.name=NULL, check.inv=FALSE, method=c('coverage', 'bipartite')){

  ## to retrieve the first sample, use something like "" in readSVvcf (NULL means all variants)
  if(is.null(sample.name)){
    sample.name = ''
  }
  ## Read/keep no-calls variants from the call set
  if(is.character(calls.gr) & length(calls.gr)==1){
    calls.gr = readSVvcf(calls.gr, keep.ins.seq=ins.seq.comp,
                         keep.ref.seq=range.seq.comp,
                         sample.name=sample.name,
                         check.inv=check.inv, nocalls=TRUE)
  }
  calls.gr = calls.gr[which(calls.gr$ac == -1)]
  if(length(calls.gr) == 0){
    return(data.frame())
  } else {
    calls.gr$ac = 1
  }
  ## Truth set
  if(is.character(truth.gr) & length(truth.gr)==1){
    truth.gr = readSVvcf(truth.gr, keep.ins.seq=ins.seq.comp,
                         keep.ref.seq=range.seq.comp,
                         sample.name=sample.name,
                         check.inv=check.inv, keep.ids=TRUE)
  }
  if(length(truth.gr) == 0){
    stop("Truth set has no SVs.")
  }
  
  ## Prepare overlaps between no-calls and truth set
  ol.gr = prepareOl(truth.gr, calls.gr, min.rol=min.del.rol,
                   max.ins.dist=max.ins.dist,
                   range.seq.comp=range.seq.comp,
                   ins.seq.comp=ins.seq.comp, nb.cores=nb.cores)

  ## Annotation overlaps
  ol.gr = annotateOl(ol.gr, min.ol=min.ol, method=method)
  
  ## Extract no-calls variants
  eval.o = evalOl(ol.gr, truth.gr, calls.gr)
  nocalls.df = lapply(names(eval.o$regions), function(svtype){
    regs = eval.o$regions[[svtype]]
    tps = regs$TP
    if(length(tps)==0){
      return(NULL)
    }
    tps = GenomicRanges::resize(tps, fix='center',
                                width=GenomicRanges::width(tps) + 2*max.ins.dist)
    tps = as.data.frame(tps)
    cols.tokeep = c('seqnames','start','end','type')
    if(!is.null(tps$svid)) cols.tokeep = c(cols.tokeep, 'svid')
    tps[, cols.tokeep]
  })
  nocalls.df = do.call(rbind, nocalls.df)
  
  return(nocalls.df)
}
