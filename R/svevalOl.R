##' @title SV evaluation based on overlap and variant size
##' @param calls.gr call set. A GRanges or the path to a VCF file.
##' @param truth.gr truth set. A GRanges or the path to a VCF file.
##' @param max.ins.dist maximum distance for insertions to be clustered. Default is 20.
##' @param min.cov the minimum coverage to be considered a match. Default is 0.5
##' @param min.del.rol minimum reciprocal overlap for deletions. Default is 0.1
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @param min.size the minimum SV size to be considered. Default 0.
##' @param max.size the maximum SV size to be considered. Default is Inf.
##' @param bed.regions If non-NULL, a GRanges object or path to a BED file
##' (no headers) with regions of interest.
##' @param bed.regions.ol minimum proportion of sv.gr that must overlap
##' regions.gr. Default is 0.5
##' @param sample.name the name of the sample to use if VCF files given as
##' input. If NULL (default), use first sample.
##' @param outfile the TSV file to output the results. If NULL (default), returns a data.frame.
##' @param out.bed.prefix prefix for the output BED files. If NULL (default), no BED output.
##' @return a data.frame with TP, FP and FN for each SV type.
##' @author Jean Monlong
##' @export
svevalOl <- function(calls.gr, truth.gr, max.ins.dist=20, min.cov=.5,
                     min.del.rol=.1, ins.seq.comp=FALSE, nb.cores=1,
                     min.size=0, max.size=Inf, bed.regions=NULL,
                     bed.regions.ol=.5,
                     sample.name=NULL, outfile=NULL, out.bed.prefix=NULL){
  if(is.character(calls.gr) & length(calls.gr)==1){
    message('Importing ', calls.gr)
    calls.gr = readSVvcf(calls.gr, keep.ins.seq=ins.seq.comp, sample.name=sample.name)
  }
  if(is.character(truth.gr) & length(truth.gr)==1){
    message('Importing ', truth.gr)
    truth.gr = readSVvcf(truth.gr, keep.ins.seq=ins.seq.comp, sample.name=sample.name)
  }
  if(min.size>0 | !is.infinite(max.size)){
    message('Filtering SVs by size.')
    calls.gr = calls.gr[which(calls.gr$size>=min.size &
                              calls.gr$size<=max.size)]
    truth.gr = truth.gr[which(truth.gr$size>=min.size &
                              truth.gr$size<=max.size)]
  }
  if(!is.null(bed.regions)){
    message('Keeping SVs overlapping regions of interest')
    if(is.character(bed.regions) & length(bed.regions) == 1){
      bed.regions = utils::read.table(bed.regions, sep='\t', as.is=TRUE)
      colnames(bed.regions)[1:3] = c('chr','start','end')
      bed.regions = GenomicRanges::makeGRangesFromDataFrame(bed.regions)
    }
    calls.gr = filterSVs(calls.gr, bed.regions, ol.prop=bed.regions.ol)
    truth.gr = filterSVs(truth.gr, bed.regions, ol.prop=bed.regions.ol)
  }
  ol.l = suppressWarnings(
    annotateOl(calls.gr, truth.gr, max.ins.gap=max.ins.dist,
               min.del.rol=min.del.rol, ins.seq.comp=ins.seq.comp,
               nb.cores=nb.cores)
  )
  eval.df = evalOl(ol.l, min.cov=min.cov, outprefix=out.bed.prefix)
  if(!is.null(outfile)){
    utils::write.table(eval.df, file=outfile, sep='\t', row.names=FALSE, quote=FALSE)
    return(outfile)
  } else {
    return(eval.df)
  }
}
