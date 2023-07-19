##' Compares SVs from a call-set to SVs from a truth-set.
##'
##' Different overlapping approaches are available. See \code{svOverlap} for more details.
##' We recommend using the default \code{method='coverage'} when evaluating the calls (absence/presence
##' of SVs) and \code{method='bipartite'} when evaluating genotypes. The latter will match
##' a variant in the call-set with at most one variant in the truth-set which will penalyze.
##' one-to-many configurations as is usually preferred when comparing exact genotypes. To further
##' switch on the 'genotype evaluation' mode, use \code{geno.eval=TRUE}. It can also help to stitch
##' fragmented calls and merge heterozygotes into homozygotes using \code{merge.hets=TRUE} and
##' \code{stitch.hets=TRUE} (see example below).
##'
##' The evaluation will be performed for different quality thresholds on the call-set in order to make
##' a precision-recall curve. If the VCF are to be read it will look for the field specified in
##' \code{qual.field} (GQ by default, and QUAL if not found). If you don't want the PR curve,
##' the evaluation can be sped up by running with for example \code{qual.ths=0}.
##' 
##' Equivalent SVs are sometimes recorded as quite different variants because placed at
##' different locations of a short tandem repeat. For example, imagine a large 100 bp
##' tandem repeat in the reference genome. An expansion of 50 bp might be represented
##' as a 50 bp insertion at the beginning of the repeat in the callset but at the end
##' of the repeat in the truth set. Because they are distant by 100 bp they might not
##' match. Instead of increasing the distance threshold too much, passing an annotation of
##' known simple repeats in the \code{simprep=} parameter provides
##' a more flexible way of matching variants by first extending them with nearby simple
##' repeats. In this example, because we know of this tandem repeat, both insertions will
##' be extended to span the full annotated reference repeat, hence ensuring that they are
##' matched and compared (e.g. by reciprocal size or sequence alignment distance)
##' short tandem repeat. 
##' @title SV evaluation based on overlap and variant size
##' @param calls.gr call set. A GRanges or the path to a VCF file.
##' @param truth.gr truth set. A GRanges or the path to a VCF file.
##' @param max.ins.dist maximum distance for insertions to be clustered. Default is 20.
##' @param min.ol the minimum overlap/coverage to be considered a match. Default is 0.5
##' @param min.del.rol minimum reciprocal overlap for deletions. Default is 0.1
##' @param range.seq.comp compare sequence instead of only overlapping deletions/inversion/etc. Default is FALSE.
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @param min.size the minimum SV size to be considered. Default 50.
##' @param max.size the maximum SV size to be considered. Default is Inf.
##' @param bed.regions If non-NULL, a GRanges object or path to a BED file
##' (no headers) with regions of interest.
##' @param bed.regions.ol minimum proportion of sv.gr that must overlap
##' regions.gr. Default is 0.5
##' @param qual.field fields to use as quality. 
##' @param sample.name the name of the sample to use if VCF files given as
##' input. If NULL (default), use first sample.
##' @param outfile the TSV file to output the results. If NULL (default), returns a data.frame.
##' @param out.bed.prefix prefix for the output BED files. If NULL (default), no BED output.
##' @param qual.ths the quality thresholds for the PR curve. If NULL, will use quantiles. see the \code{qual.quantiles} parameter below.
##' @param qual.quantiles the quality quantiles for the PR curve, if qual.ths is NULL. Default is (0, .1, ..., .9, 1).
##' @param check.inv should the sequence of MNV be compared to identify inversions. 
##' @param geno.eval should het/hom be evaluated separately (genotype evaluation). Default
##' FALSE.
##' @param stitch.hets should clustered hets be stitched together before genotype evatuation.
##' Default is FALSE.
##' @param stitch.dist the maximum distance to stitch hets during genotype evaluation.
##' @param merge.hets should similar hets be merged into homs before genotype evaluation.
##' Default is FALSE.
##' @param merge.rol the minimum reciprocal overlap to merge hets before genotype
##' evaluation.
##' @param method the method to annotate the overlap. Either 'coverage' (default) for the
##' cumulative coverage (e.g. to deal with fragmented calls); or 'bipartite' for a 1-to-1
##' matching of variants in the calls and truth sets.
##' @param simprep optional simple repeat annotation. Default is NULL. If non-NULL, GRanges to be used to
##' extend variants when overlapping/clustering
##' @param log.level the level of information in the log. Default is "CRITICAL" (basically no log).
##' @return a list with
##' \item{eval}{a data.frame with TP, FP and FN for each SV type when including all variants}
##' \item{curve}{a data.frame with TP, FP and FN for each SV type when using different quality thesholds}
##' \item{svs}{a list of GRanges object with FP, TP and FN for each SV type (using quality threshold with best F1).}
##' \item{mqual.bestf1}{the quality threshold that produces best F1 (and corresponding to 'svs' GRanges).}
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
##'
##' ## Genotype evaluation
##' eval = svevalOl(calls.gr, truth.gr, geno.eval=TRUE, merge.hets=TRUE,
##'                 stitch.hets=TRUE, method='bipartite')
##' }
svevalOl <- function(calls.gr, truth.gr, max.ins.dist=20, min.ol=.5,
                     min.del.rol=.1, range.seq.comp=FALSE, ins.seq.comp=FALSE, nb.cores=1,
                     min.size=50, max.size=Inf, bed.regions=NULL,
                     bed.regions.ol=.5, qual.field=c('GQ', 'QUAL'),
                     sample.name=NULL, outfile=NULL,
                     out.bed.prefix=NULL,
                     qual.ths=c(0, 2, 3, 4, 5, 7, 10, 12, 14, 21, 27, 35, 45, 50, 60, 75, 90, 99, 110, 133, 167, 180, 250, 350, 450, 550, 650),
                     qual.quantiles=seq(0,1,.1),
                     check.inv=FALSE, geno.eval=FALSE, stitch.hets=FALSE,
                     stitch.dist=20, merge.hets=FALSE, merge.rol=.8, method=c('coverage', 'bipartite'),
                     simprep=NULL,
                     log.level=c('CRITICAL', 'WARNING', 'INFO')){
  logging::setLevel(log.level[1])
  ## to retrieve the first sample, use something like "" in readSVvcf (NULL means all variants)
  if(is.null(sample.name)){
    sample.name = ''
  }
  if(is.character(calls.gr) & length(calls.gr)==1){
    logging::loginfo(paste('Reading SVs from', calls.gr))
    calls.gr = readSVvcf(calls.gr, keep.ref.seq=range.seq.comp, keep.ins.seq=ins.seq.comp,
                         qual.field=qual.field, sample.name=sample.name, check.inv=check.inv)
  }
  if(is.character(truth.gr) & length(truth.gr)==1){
    logging::loginfo(paste('Reading SVs from', truth.gr))
    truth.gr = readSVvcf(truth.gr, keep.ref.seq=range.seq.comp, keep.ins.seq=ins.seq.comp,
                         qual.field=qual.field, sample.name=sample.name, check.inv=check.inv)
  }
  # keep only non-ref calls
  calls.gr = calls.gr[which(calls.gr$ac>0)]
  truth.gr = truth.gr[which(truth.gr$ac>0)]
  if(length(truth.gr) == 0){
    stop("Truth set has no SVs.")
  }
  if(length(calls.gr) == 0){
    warning("Call set has no SVs.")
  }
  if(!is.null(bed.regions)){
    if(is.character(bed.regions) & length(bed.regions) == 1){
      logging::loginfo(paste('Reading regions from', bed.regions))
      bed.regions = utils::read.table(bed.regions, sep='\t', as.is=TRUE)
      colnames(bed.regions)[1:3] = c('chr','start','end')
      bed.regions = GenomicRanges::makeGRangesFromDataFrame(bed.regions)
    }
  }

  ## if BND and TRA, homogeneize the SV type
  if(any(calls.gr$type %in% c('BND', 'TRA'))){
    calls.gr$type = ifelse(calls.gr$type %in% c('BND', 'TRA'), 'BND', calls.gr$type)
  }
  if(any(truth.gr$type %in% c('BND', 'TRA'))){
    truth.gr$type = ifelse(truth.gr$type %in% c('BND', 'TRA'), 'BND', truth.gr$type)
  }
  ## optional: extend variants with simple repeat annotation
  if(!is.null(simprep)){
    logging::loginfo('Extend callset using simple repeat annotation')
    calls.gr = extendSVwithSimpRep(calls.gr, simprep)
    logging::loginfo('Extend truthset using simple repeat annotation')
    truth.gr = extendSVwithSimpRep(truth.gr, simprep)
  }
  
  ## If evaluation per genotype, do we want to stitch and/or merge heterozygous variants?
  if(length(calls.gr)>0 & length(truth.gr)>0){
    if(geno.eval & (stitch.hets | merge.hets)){
      calls.gr = stitchMergeHets(calls.gr, do.stitch=stitch.hets, do.merge=merge.hets, min.rol=merge.rol,
                                 max.ins.dist=max.ins.dist,
                                 range.seq.comp=range.seq.comp, ins.seq.comp=ins.seq.comp, stitch.dist=stitch.dist)
      truth.gr = stitchMergeHets(truth.gr, do.stitch=stitch.hets, do.merge=merge.hets, min.rol=merge.rol,
                                 max.ins.dist=max.ins.dist,
                                 range.seq.comp=range.seq.comp, ins.seq.comp=ins.seq.comp,
                                 stitch.dist=stitch.dist)
    } 
  }

  ## Prepare the overlaps, by SV type and eventually by genotype
  ol.gr = prepareOl(truth.gr, calls.gr, min.rol=min.del.rol,
                   max.ins.dist=max.ins.dist,
                   range.seq.comp=range.seq.comp, ins.seq.comp=ins.seq.comp,
                   nb.cores=nb.cores,
                   by.gt=geno.eval)

  ## Filter SVs out of the desired size range or regions
  calls.gr = filterSVs(calls.gr, regions.gr=bed.regions,
                       ol.prop=bed.regions.ol,
                       min.size=min.size, max.size=max.size,
                       mark.pass=TRUE)
  truth.gr = filterSVs(truth.gr, regions.gr=bed.regions,
                       ol.prop=bed.regions.ol,
                       min.size=min.size, max.size=max.size,
                       mark.pass=TRUE)
  if(length(ol.gr)>0){
    ol.gr = ol.gr[which((ol.gr$queryHits %in% which(truth.gr$pass)) |
                        (ol.gr$subjectHits %in% which(calls.gr$pass)))]
  }

  ## Compute evaluation metrics for each quality threshold
  if(!is.null(qual.ths)){
    qual.r = unique(c(0, qual.ths))
  } else {
    qual.r = unique(c(0, stats::quantile(calls.gr$qual, probs=qual.quantiles)))
  }
  eval.quals.o = parallel::mclapply(qual.r, function(mqual){
    ## message(mqual)
    ## remove calls with quality lower than threshold
    calls.gr$pass[which(calls.gr$qual < mqual)] = FALSE
    ol.gr = ol.gr[which(ol.gr$subjectHits %in% which(calls.gr$pass))]
    ## annotate the overlaps and compute evaluation metrics
    ol.gr = annotateOl(ol.gr, min.ol=min.ol, method=method)
    eval.o = evalOl(ol.gr, truth.gr, calls.gr)
    eval.o$eval$qual = mqual
    eval.o
  }, mc.cores=nb.cores)
  ## all metrics in a data.frame
  eval.curve.df = do.call(rbind, lapply(eval.quals.o, function(ll) ll$eval))
  ## best quality threshold (maximizing total F1)
  f1s = sapply(eval.quals.o, function(ll) ll$eval$F1[which(ll$eval$type=='Total')])
  bestf1 = utils::head(order(f1s, decreasing=TRUE), 1)
  eval.bestf1 = eval.quals.o[[bestf1]]
  mqual.bestf1 = qual.r[bestf1]
    
  ## Write BED files with FP, TP, FN
  if(!is.null(out.bed.prefix)){
    tmp = lapply(names(eval.bestf1$regions), function(svtype){
      regs = eval.bestf1$regions[[svtype]]
      utils::write.table(as.data.frame(regs$FN), file=paste0(out.bed.prefix, svtype, '-FN.tsv'), sep='\t', row.names=TRUE, quote=FALSE)
      utils::write.table(as.data.frame(regs$TP.baseline), file=paste0(out.bed.prefix, svtype, '-TP-baseline.tsv'), sep='\t', row.names=TRUE, quote=FALSE)
      utils::write.table(as.data.frame(regs$FP), file=paste0(out.bed.prefix, svtype, '-FP.tsv'), sep='\t', row.names=TRUE, quote=FALSE)
      utils::write.table(as.data.frame(regs$TP), file=paste0(out.bed.prefix, svtype, '-TP-call.tsv'), sep='\t', row.names=TRUE, quote=FALSE)
    })
  }

  ## Write results for PR curve
  if(!is.null(out.bed.prefix)){
    utils::write.table(eval.curve.df, file=paste0(out.bed.prefix, 'prcurve.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
    grDevices::pdf(paste0(out.bed.prefix, 'prcurve.pdf'), 9, 7)
    tmp = lapply(plot_prcurve(paste0(out.bed.prefix, 'prcurve.tsv')), print)
    grDevices::dev.off()
  }
  
  if(!is.null(outfile)){
    utils::write.table(eval.bestf1$eval, file=outfile, sep='\t', row.names=FALSE, quote=FALSE)
  }
  return(list(eval=eval.bestf1$eval, curve=eval.curve.df, svs=eval.bestf1$regions, mqual.bestf1=mqual.bestf1))
}
