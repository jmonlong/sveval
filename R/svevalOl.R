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
##' @param qual.quantiles the QUAL quantiles for the PR curve. Default is (0, .1, ..., .9, 1).
##' @return a list with
##' \item{eval}{a data.frame with TP, FP and FN for each SV type when including all variants}
##' \item{curve}{a data.frame with TP, FP and FN for each SV type when using different quality thesholds}
##' @author Jean Monlong
##' @export
##' @examples
##' \dontrun{
##' ## From VCF files
##' eval = svevalOl('calls.vcf', 'truth.vcf')
##'
##' ## From GRanges
##' calls.gr = readSVvcf('calls.vcf')
##' truth.gr = readSVvcf('truth.vcf')
##' eval = svevalOl(calls.gr, truth.gr)
##' }
svevalOl <- function(calls.gr, truth.gr, max.ins.dist=20, min.cov=.5,
                     min.del.rol=.1, ins.seq.comp=FALSE, nb.cores=1,
                     min.size=0, max.size=Inf, bed.regions=NULL,
                     bed.regions.ol=.5, sample.name=NULL, outfile=NULL,
                     out.bed.prefix=NULL, qual.quantiles=seq(0,1,.1)){
  if(is.character(calls.gr) & length(calls.gr)==1){
    calls.gr = readSVvcf(calls.gr, keep.ins.seq=ins.seq.comp, sample.name=sample.name)
  }
  if(is.character(truth.gr) & length(truth.gr)==1){
    truth.gr = readSVvcf(truth.gr, keep.ins.seq=ins.seq.comp, sample.name=sample.name)
  }
  if(length(calls.gr)>0 & length(truth.gr)>0 & !is.null(bed.regions)){
    if(is.character(bed.regions) & length(bed.regions) == 1){
      bed.regions = utils::read.table(bed.regions, sep='\t', as.is=TRUE)
      colnames(bed.regions)[1:3] = c('chr','start','end')
      bed.regions = GenomicRanges::makeGRangesFromDataFrame(bed.regions)
    }
  }

  ## Overlap SVs
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

  ## Compute coverage and evaluation metrics
  qual.r = unique(c(0, stats::quantile(calls.gr$QUAL, probs=qual.quantiles)))
  eval.curve.df = lapply(qual.r, function(mqual){
    ins.a = annotateOl(ol.ins, min.qual=mqual)
    del.a = annotateOl(ol.del, min.qual=mqual)
    inv.a = annotateOl(ol.inv, min.qual=mqual)
    ol.l = list(
      calls=c(ins.a$calls, del.a$calls, inv.a$calls),
      truth=c(ins.a$truth, del.a$truth, inv.a$truth)
    )
    if(length(ol.l$calls)==0 | length(ol.l$truth)==0){
      eval.df = evalOl(NULL)
    } else {
      ol.l$calls = filterSVs(ol.l$calls, regions.gr=bed.regions, ol.prop=bed.regions.ol,
                             min.size=min.size, max.size=max.size)
      ol.l$truth = filterSVs(ol.l$truth, regions.gr=bed.regions, ol.prop=bed.regions.ol,
                             min.size=min.size, max.size=max.size)
      op = out.bed.prefix
      if(mqual>0){
        op=NULL
      }
      eval.df = evalOl(ol.l, min.cov=min.cov, outprefix=op)
    }
    eval.df$qual = mqual
    eval.df
  })
  eval.curve.df = do.call(rbind, eval.curve.df)
  eval.df = eval.curve.df[which(eval.curve.df$qual==0),]
  eval.df$qual = NULL

  ## Write results for PR curve
  if(!is.null(out.bed.prefix)){
    utils::write.table(eval.curve.df, file=paste0(out.bed.prefix, 'prcurve.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
    grDevices::pdf(paste0(out.bed.prefix, 'prcurve.pdf'), 9, 7)
    tmp = lapply(plot_prcurve(paste0(out.bed.prefix, 'prcurve.tsv')), print)
    grDevices::dev.off()
  }
  
  if(!is.null(outfile)){
    utils::write.table(eval.df, file=outfile, sep='\t', row.names=FALSE, quote=FALSE)
    return(outfile)
  } else {
    return(list(eval=eval.df, curve=eval.curve.df))
  }
}
