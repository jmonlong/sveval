##' @title SV evaluation based on overlap and variant size
##' @param calls.gr call set. A GRanges or the path to a VCF file.
##' @param truth.gr truth set. A GRanges or the path to a VCF file.
##' @param max.ins.dist maximum distance for insertions to be clustered. Default is 20.
##' @param min.cov the minimum coverage to be considered a match. Default is 0.5
##' @param min.del.rol minimum reciprocal overlap for deletions. Default is 0.1
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param nb.cores number of processors to use. Default is 1.
##' @param min.size the minimum SV size to be considered. Default 50.
##' @param max.size the maximum SV size to be considered. Default is Inf.
##' @param bed.regions If non-NULL, a GRanges object or path to a BED file
##' (no headers) with regions of interest.
##' @param bed.regions.ol minimum proportion of sv.gr that must overlap
##' regions.gr. Default is 0.5
##' @param qual.field fields to use as quality. Will be tried in order.
##' @param sample.name the name of the sample to use if VCF files given as
##' input. If NULL (default), use first sample.
##' @param outfile the TSV file to output the results. If NULL (default), returns a data.frame.
##' @param out.bed.prefix prefix for the output BED files. If NULL (default), no BED output.
##' @param qual.ths the QUAL thresholds for the PR curve. If NULL, will use quantiles. see \code{qual.quantiles}.
##' @param qual.quantiles the QUAL quantiles for the PR curve, if qual.ths is NULL. Default is (0, .1, ..., .9, 1).
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
##'
##' ## Genotype evaluation
##' eval = svevalOl(calls.gr, truth.gr, geno.eval=TRUE, merge.hets=TRUE, stitch.hets=TRUE)
##' }
svevalOl <- function(calls.gr, truth.gr, max.ins.dist=20, min.cov=.5,
                     min.del.rol=.1, ins.seq.comp=FALSE, nb.cores=1,
                     min.size=50, max.size=Inf, bed.regions=NULL,
                     bed.regions.ol=.5, qual.field=c('QUAL', 'GQ'),
                     sample.name=NULL, outfile=NULL,
                     out.bed.prefix=NULL,
                     qual.ths=c(0, 2, 3, 4, 5, 7, 10, 12, 14, 21, 27, 35, 45, 50, 60, 75, 90, 99, 110, 133, 167, 180, 250, 350, 450, 550, 650),
                     qual.quantiles=seq(0,1,.1),
                     check.inv=FALSE, geno.eval=FALSE, stitch.hets=FALSE,
                     stitch.dist=20, merge.hets=FALSE, merge.rol=.8){
  if(is.character(calls.gr) & length(calls.gr)==1){
    calls.gr = readSVvcf(calls.gr, keep.ins.seq=ins.seq.comp, qual.field=qual.field, sample.name=sample.name, check.inv=check.inv)
  }
  if(is.character(truth.gr) & length(truth.gr)==1){
    truth.gr = readSVvcf(truth.gr, keep.ins.seq=ins.seq.comp, qual.field=qual.field, sample.name=sample.name, check.inv=check.inv)
  }
  if(length(truth.gr) == 0){
    stop("Truth set has no SVs.")
  }
  if(length(calls.gr)>0 & length(truth.gr)>0 & !is.null(bed.regions)){
    if(is.character(bed.regions) & length(bed.regions) == 1){
      bed.regions = utils::read.table(bed.regions, sep='\t', as.is=TRUE)
      colnames(bed.regions)[1:3] = c('chr','start','end')
      bed.regions = GenomicRanges::makeGRangesFromDataFrame(bed.regions)
    }
  }

  ## If not per genotype, set every variant to homozygous
  if(length(calls.gr)>0 & length(truth.gr)>0){
    if(geno.eval){
      iterStitch <- function(svs.gr, stitch.dist){
        svs.gr = lapply(unique(svs.gr$type), function(type){
          svs.t = svs.gr[which(svs.gr$type==type)]
          nhets = Inf
          while(length((hets.idx = which(svs.t$GT == 'het'))) < nhets){
            nhets = length(hets.idx)
            hets = stitchSVs(svs.t[hets.idx], stitch.dist=stitch.dist)
            svs.t = c(hets, svs.t[which(svs.t$GT == 'hom')])
          }
          return(svs.t)
        })
        do.call(c, svs.gr)
      }
      iterMerge <- function(svs.gr, min.rol, max.ins.gap, ins.seq.comp){
        svs.gr = lapply(unique(svs.gr$type), function(type){
          svs.t = svs.gr[which(svs.gr$type==type)]
          nhets = Inf
          while(length((hets.idx = which(svs.t$GT == 'het'))) < nhets){
            nhets = length(hets.idx)
            hets = mergeHets(svs.t[hets.idx], min.rol=min.rol,
                             max.ins.gap=max.ins.gap, ins.seq.comp=ins.seq.comp)
            svs.t = c(hets, svs.t[which(svs.t$GT == 'hom')])
          }
          return(svs.t)
        })
        do.call(c, svs.gr)
      }
      ## Stitch hets SVs
      if(stitch.hets){
        ## Merge hets once first
        if(merge.hets){
          calls.gr = iterMerge(calls.gr, min.rol=merge.rol,
                               max.ins.gap=max.ins.dist,
                               ins.seq.comp=ins.seq.comp)
          truth.gr = iterMerge(truth.gr, min.rol=merge.rol,
                               max.ins.gap=max.ins.dist,
                               ins.seq.comp=ins.seq.comp)
        }
        calls.gr = iterStitch(calls.gr, stitch.dist=stitch.dist)
        truth.gr = iterStitch(truth.gr, stitch.dist=stitch.dist)
      }
      ## Merge hets
      if(merge.hets){
        calls.gr = iterMerge(calls.gr, min.rol=merge.rol,
                             max.ins.gap=max.ins.dist,
                             ins.seq.comp=ins.seq.comp)
        truth.gr = iterMerge(truth.gr, min.rol=merge.rol,
                             max.ins.gap=max.ins.dist,
                             ins.seq.comp=ins.seq.comp)
      }
    } else {
      truth.gr$GT = 'hom'
      calls.gr$GT = 'hom'
    } 
  }
  
  ## Overlap per genotype
  ol.gt = lapply(unique(c(truth.gr$GT, calls.gr$GT)), function(gt){
    calls.gr = calls.gr[which(calls.gr$GT == gt)]
    truth.gr = truth.gr[which(truth.gr$GT == gt)]
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
    return(list(ol.ins=ol.ins, ol.del=ol.del, ol.inv=ol.inv))
  })

  ## Compute coverage and evaluation metrics
  if(!is.null(qual.ths)){
    qual.r = unique(c(0, qual.ths))
  } else {
    qual.r = unique(c(0, stats::quantile(calls.gr$QUAL, probs=qual.quantiles)))
  }
  eval.curve.df = lapply(qual.r, function(mqual){
    ## Insertion annotation for each genotype
    ins.a.gt = lapply(ol.gt, function(ll) annotateOl(ll$ol.ins, min.qual=mqual))
    ## Deletion annotation for each genotype
    del.a.gt = lapply(ol.gt, function(ll) annotateOl(ll$ol.del, min.qual=mqual))
    ## Inversion annotation for each genotype
    inv.a.gt = lapply(ol.gt, function(ll) annotateOl(ll$ol.inv, min.qual=mqual))

    ol.l = c(ins.a.gt, del.a.gt, inv.a.gt)
    ol.l = list(
      calls=do.call(c, lapply(ol.l, function(ll) ll$calls)),
      truth=do.call(c, lapply(ol.l, function(ll) ll$truth))
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
