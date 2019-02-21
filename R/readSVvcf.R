##' Read a VCF file that contains SVs and create a GRanges with relevant information, e.g. SV size or genotype quality.
##'
##' By default, the quality information is taken from the QUAL field. If all
##' values are NA or 0, the function will try other fields as speficied in
##' the "qual.field" vector. Fields can be from the INFO or FORMAT fields.
##' @title Read SVs from a VCF file
##' @param vcf.file the path to the VCF file
##' @param keep.ins.seq should it keep the inserted sequence? Default is FALSE.
##' @param sample.name the name of the sample to use. If NULL (default), use
##' first sample.
##' @param qual.field fields to use as quality. Will be tried in order.
##' @return a GRanges object with relevant information.
##' @author Jean Monlong
##' @export
##' @examples
##' \dontrun{
##' calls.gr = readSVvcf('calls.vcf')
##' }
readSVvcf <- function(vcf.file, keep.ins.seq=FALSE, sample.name=NULL, qual.field=c('QUAL', 'GQ')){
  vcf = VariantAnnotation::readVcf(vcf.file, row.names=FALSE)
  gr = DelayedArray::rowRanges(vcf)
  ## If sample specified, retrieve appropriate GT
  GT.idx = 1
  if(!is.null(sample.name)){
    GT.idx = which(sample.name == colnames(VariantAnnotation::geno(vcf)$GT))
  } 
  gr$GT = unlist(VariantAnnotation::geno(vcf)$GT[, GT.idx])
  
  ## Remove obvious SNVs
  singlealt = which(unlist(lapply(gr$ALT, length))==1)
  alt.sa = unlist(Biostrings::nchar(gr$ALT[singlealt]))
  ref.sa = Biostrings::nchar(gr$REF[singlealt])
  snv.idx = singlealt[which(ref.sa==1 & alt.sa==1)]
  snv.idx = snv.idx[which(as.character(unlist(gr$ALT[snv.idx])) != as.character(gr$REF[snv.idx]))]
  if(length(snv.idx)>0){
    nonsnv.idx = setdiff(1:length(gr), snv.idx)
    gr = gr[nonsnv.idx]
    vcf = vcf[nonsnv.idx]
  }
  
  ## Convert missing qualities to 0
  if(any(is.na(gr$QUAL))){
    gr$QUAL[which(is.na(gr$QUAL))] = 0
  }
  ## Extract quality information
  qual.found = FALSE
  qfield.ii = 1
  while(!qual.found & qfield.ii < length(qual.field) + 1){
    if(qual.field[qfield.ii] == 'QUAL' & any(gr$QUAL>0)){
      qual.found = TRUE
    } else if(qual.field[qfield.ii] %in% names(VariantAnnotation::geno(vcf))){
      gr$QUAL = unlist(VariantAnnotation::geno(vcf)[[qual.field[qfield.ii]]][, GT.idx])
      qual.found = TRUE
    } else if(qual.field[qfield.ii] %in% colnames(VariantAnnotation::info(vcf))){
      gr$QUAL = unlist(VariantAnnotation::info(vcf)[[qual.field[qfield.ii]]])
      qual.found = TRUE
    }
    qfield.ii = qfield.ii + 1
  }

  ## Symbolic alleles or ALT/REF ?
  if(all(c('SVTYPE', 'SVLEN') %in% colnames(VariantAnnotation::info(vcf)))){
    ## Symbolic alleles
    gr$type = unlist(VariantAnnotation::info(vcf)$SVTYPE)
    gr$size = abs(unlist(VariantAnnotation::info(vcf)$SVLEN))
    GenomicRanges::end(gr) = ifelse(gr$type == 'DEL',
                                    GenomicRanges::end(gr) + gr$size,
                                    GenomicRanges::end(gr))
  } else {
    ## ALT/REF
    ## Define variant type from alt/ref size
    alt.s = unlist(lapply(Biostrings::nchar(gr$ALT), max))
    ## Using largest ALT allele !
    ## Maybe we can do better.
    ref.s = Biostrings::nchar(gr$REF)
    gr$type = ifelse(alt.s>ref.s, 'INS', 'DEL')
    gr$type = ifelse(alt.s==ref.s, 'MNV', gr$type)
    gr$type = ifelse(alt.s==1 & ref.s==1, 'SNV', gr$type)
    ## Variants other than clear DEL, INS or SNV. 
    others = which(alt.s>1 & ref.s>1)
    if(length(others)>0){
      gr.inv = gr[others]
      alt.seq = lapply(gr.inv$ALT, function(alt)alt[which.max(Biostrings::nchar(alt))])
      alt.seq = do.call(c, alt.seq)
      ref.seq = gr.inv$REF
      isinv = checkInvSeq(ref.seq, alt.seq)
      gr$type[others] = ifelse(isinv, 'INV', gr$type[others])      
    }    
    gr$size = ifelse(gr$type=='INS', alt.s, GenomicRanges::width(gr))
  }
  ## read support if available
  if('AD' %in% rownames(VariantAnnotation::geno(VariantAnnotation::header(vcf)))){
    ad.l = VariantAnnotation::geno(vcf)$AD[, GT.idx]
    gr$ref.cov = unlist(lapply(ad.l, '[', 1))
    gr$alt.cov = unlist(lapply(ad.l, '[', 2))
  } else if(all(c('RO', 'AO') %in% rownames(VariantAnnotation::geno(VariantAnnotation::header(vcf))))){
    gr$ref.cov = VariantAnnotation::geno(vcf)$RO
    gr$alt.cov = unlist(lapply(VariantAnnotation::geno(vcf)$AO, '[', 1))
  } else {
    gr$alt.cov = gr$ref.cov = NA
  }
  ## Remove unused columns
  gr$REF = gr$paramRangeID = gr$FILTER = NULL
  if(!keep.ins.seq){
    gr$ALT = NULL
  }
  ## Remove "ref" variants and SNVs
  gr = gr[which(gr$GT!='0' & gr$GT!='0/0'  & gr$GT!='0|0' &
                gr$GT!='./.' & gr$GT!='.')]
  gr = gr[which(gr$type!='SNV' & gr$type!='MNV')]
  return(gr)
}
