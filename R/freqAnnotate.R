##' Input SVs are matched to SVs in the catalog. Each SV is then annotated with the
##' maximum frequency across all the variants in the catalog that match. Although the
##' different overlap approaches could be used here, the 'reciprocal' overlap makes the
##' more sense because we want to list all possible matches without penalizing over-matching
##' (like 'bipartite' would) or include fragmented calls (like 'coverage' would) as they
##' could be very variants with uncomparable frequencies.
##'
##' @title Annotate SVs with frequency in catalog
##' @param svs a VCF object with SVs to annotate.
##' @param cat a VCF object with the SV catalog with frequency estimates.
##' @param min.ol the minimum overlap/coverage to be considered a match. Default is 0.5
##' @param min.del.rol minimum reciprocal overlap for deletions. Default is 0.1
##' @param max.ins.dist maximum distance for insertions to be clustered. Default is 20.
##' @param check.inv should the sequence of MNV be compared to identify inversions. 
##' @param range.seq.comp compare sequence instead of only overlapping deletions/inversions/etc. Default is FALSE.
##' @param ins.seq.comp compare sequence instead of insertion sizes. Default is FALSE.
##' @param out.vcf If non-NULL, write output to this VCF file. 
##' @param freq.field the field with the frequency estimate in the 'cat' input. Default is 'AF'.
##' @param out.freq.field the new field's name. Default is 'AFMAX'
##' @param method the method to annotate the overlap. Recommended is 'reciprocal' (default). See details.
##' @param nb.cores number of processors to use. Default is 1.
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
##' calls.vcf = readSVvcf('calls.vcf', out.fmt="vcf")
##' cat.vcf = readSVvcf('gnomad.vcf', out.fmt="vcf")
##' calls.freq.vcf = freqAnnotate(calls.vcf, cat.vcf)
##' }
freqAnnotate <- function(svs, cat, min.ol=.5, min.del.rol=.1, max.ins.dist=20, check.inv=FALSE,
                         range.seq.comp=FALSE, ins.seq.comp=FALSE, out.vcf=NULL, freq.field='AF',
                         out.freq.field='AFMAX', method=c('reciprocal', 'coverage', 'bipartite'),
                         nb.cores=1){

  if(is.character(svs) && length(svs) == 1){
    svs = readSVvcf(svs, out.fmt='vcf', check.inv=check.inv)
  }
  if(is.character(cat) && length(cat) == 1){
    cat = readSVvcf(cat, out.fmt='vcf', check.inv=check.inv, other.field=freq.field)
  }
  
  ## Read/parse inputs
  svs.gr = DelayedArray::rowRanges(svs)
  svs.gr$size = VariantAnnotation::info(svs)$SVLEN
  svs.gr$type = VariantAnnotation::info(svs)$SVTYPE
  GenomicRanges::end(svs.gr) = VariantAnnotation::info(svs)$END
  cat.gr = DelayedArray::rowRanges(cat)
  cat.gr$size = VariantAnnotation::info(cat)$SVLEN
  cat.gr$type = VariantAnnotation::info(cat)$SVTYPE

  desc = VariantAnnotation::info(VariantAnnotation::header(cat))
  if(!(freq.field %in% rownames(desc))){
    stop(freq.field, ' is not an INFO field in the catalog')
  }
  if(desc[freq.field, 'Number'] == '1'){
    cat.gr$freq = VariantAnnotation::info(cat)[[freq.field]]
  } else if(desc[freq.field, 'Number'] == 'A'){
    cat.gr$freq = unlist(lapply(VariantAnnotation::info(cat)[[freq.field]], sum))
  }
  
  GenomicRanges::end(cat.gr) = VariantAnnotation::info(cat)$END

  ## Prepare VCF field
  freqs = rep(0, length(svs))

  ## Overlap SVs
  ol.gr = prepareOl(cat.gr, svs.gr, min.rol=min.del.rol, max.ins.dist=max.ins.dist,
                    range.seq.comp=range.seq.comp, ins.seq.comp=ins.seq.comp,
                    nb.cores=nb.cores, by.gt=FALSE)
  ol.gr = annotateOl(ol.gr, min.ol=min.ol, method=method)

  if(length(ol.gr)>0){
    ## maximum frequency for each input SV
    freq.df = ol.gr %>%
      as.data.frame %>% 
      dplyr::mutate(freq=cat.gr$freq[.data$queryHits]) %>% 
      dplyr::group_by(.data$subjectHits) %>%
      dplyr::summarize(freq=max(.data$freq))
    freqs[freq.df$subjectHits] = freq.df$freq
  }
  
  ## Add frequency field
  new.info.h = S4Vectors::DataFrame(Number='1', Type='Float',
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
