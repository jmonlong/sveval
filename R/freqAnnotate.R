##' @title Annotate SVs with frequency in catalog
##' @param svs a VCF object with SVs to annotate.
##' @param cat a VCF object with the SV catalog with frequency estimates.
##' @param min.cov the minimum coverage to be considered a match. Default is 0.5
##' @param min.del.rol minimum reciprocal overlap for deletions. Default is 0.1
##' @param max.ins.dist maximum distance for insertions to be clustered. Default is 20.
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param out.vcf If non-NULL, write output to this VCF file. 
##' @param freq.field the field with the frequency estimate in the 'cat' input. Default is 'AF'.
##' @param out.freq.field the new field's name. Default is 'AFMAX'
##' @return a GRanges object.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
##' @examples
##' \dontrun{
##' ## From VCF files with output written to VCF file
##' freqAnnotate('calls.vcf', 'gnomad.vcf', out.vcf='calls.withFreq.vcf')
##'
##' ## Within R
##' calls.vcf = readSVvcf('calls.vcf', vcf.object=TRUE)
##' cat.vcf = readSVvcf('gnomad.vcf', vcf.object=TRUE)
##' calls.freq.vcf = freqAnnotate(calls.vcf, cat.vcf)
##' }
freqAnnotate <- function(svs, cat, min.cov=.5, min.del.rol=.1, max.ins.dist=20,
                         ins.seq.comp=FALSE, out.vcf=NULL, freq.field='AF',
                         out.freq.field='AFMAX'){

  if(is.character(svs) && length(svs) == 1){
    svs = readSVvcf(svs, vcf.object=TRUE)
  }
  if(is.character(cat) && length(cat) == 1){
    cat = readSVvcf(cat, vcf.object=TRUE)
  }
  
  ## Read/parse inputs
  svs.gr = DelayedArray::rowRanges(svs)
  svs.gr$size = VariantAnnotation::info(svs)$SIZE
  svs.gr$type = VariantAnnotation::info(svs)$SVTYPE
  svs.gr$QUAL = VariantAnnotation::info(svs)$QUAL
  GenomicRanges::end(svs.gr) = VariantAnnotation::info(svs)$END
  svs.gr$id = 1:length(svs.gr)
  cat.gr = DelayedArray::rowRanges(cat)
  cat.gr$size = VariantAnnotation::info(cat)$SIZE
  cat.gr$type = VariantAnnotation::info(cat)$SVTYPE
  cat.gr$QUAL = VariantAnnotation::info(cat)$QUAL

  desc = VariantAnnotation::info(VariantAnnotation::header(cat))
  if(desc[freq.field, 'Number'] == '1'){
    cat.gr$freq = VariantAnnotation::info(cat)[[freq.field]]
  } else if(desc[freq.field, 'Number'] == 'A'){
    cat.gr$freq = unlist(lapply(VariantAnnotation::info(cat)[[freq.field]], sum))
  }
  
  GenomicRanges::end(cat.gr) = VariantAnnotation::info(cat)$END
  cat.gr$id = 1:length(cat.gr)

  ## Prepare VCF field
  freqs = rep(0, length(svs))
 
  ## Overlap insertions
  ol.ins = suppressWarnings(
    olInsertions(svs.gr, cat.gr, max.ins.gap=max.ins.dist,
                 ins.seq.comp=ins.seq.comp)
  )

  ## Overlap deletions
  ol.del = suppressWarnings(
    olRanges(svs.gr, cat.gr, min.rol=min.del.rol, type='DEL')
  )

  ## Overlap inversions
  ol.inv = suppressWarnings(
    olRanges(svs.gr, cat.gr, min.rol=min.del.rol, type='INV')
  )

  ## Overlap duplications
  ol.dup = suppressWarnings(
    olRanges(svs.gr, cat.gr, min.rol=min.del.rol, type='DUP')
  )

  ## Annotate SVs
  if(!is.null(ol.ins$ol)){
    freq.ins = ol.ins$ol %>%
      dplyr::mutate(freq=ol.ins$truth$freq[.data$truth.idx],
                    id=ol.ins$calls$id[.data$call.idx]) %>% 
      dplyr::group_by(.data$id) %>%
      dplyr::summarize(freq=max(.data$freq))
    freqs[freq.ins$id] = freq.ins$freq
  }
  if(!is.null(ol.del$rol.gr)){
    freq.del = ol.del$rol.gr %>%
      as.data.frame %>% 
      dplyr::mutate(freq=ol.del$truth$freq[.data$truth.idx],
                    id=ol.del$calls$id[.data$call.idx]) %>% 
      dplyr::group_by(.data$id) %>%
      dplyr::summarize(freq=max(.data$freq))
    freqs[freq.del$id] = freq.del$freq
  }
  if(!is.null(ol.inv$rol.gr)){
    freq.inv = ol.inv$rol.gr %>%
      as.data.frame %>% 
      dplyr::mutate(freq=ol.inv$truth$freq[.data$truth.idx],
                    id=ol.inv$calls$id[.data$call.idx]) %>% 
      dplyr::group_by(.data$id) %>%
      dplyr::summarize(freq=max(.data$freq))
    freqs[freq.inv$id] = freq.inv$freq
  }
  if(!is.null(ol.dup$rol.gr)){
    freq.dup = ol.dup$rol.gr %>%
      as.data.frame %>% 
      dplyr::mutate(freq=ol.dup$truth$freq[.data$truth.idx],
                    id=ol.dup$calls$id[.data$call.idx]) %>% 
      dplyr::group_by(.data$id) %>%
      dplyr::summarize(freq=max(.data$freq))
    freqs[freq.dup$id] = freq.dup$freq
  }

  ## Add frequency field
  new.info.h = S4Vectors::DataFrame(Number='1', Type='float',
                                    Description=desc[freq.field,'Description'])
  row.names(new.info.h) = out.freq.field
  VariantAnnotation::info(VariantAnnotation::header(svs)) =
    rbind(VariantAnnotation::info(VariantAnnotation::header(svs)),
          new.info.h)
  VariantAnnotation::info(svs)[[out.freq.field]] = freqs
  
  if(!is.null(out.vcf)){
    VariantAnnotation::writeVcf(svs, file=out.vcf)
    return(out.vcf)
  } else {
    return(svs)
  }
}
