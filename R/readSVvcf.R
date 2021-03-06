##' Read a VCF file that contains SVs and create a GRanges with relevant information, e.g. SV size or genotype quality.
##' 
##' By default, the quality information is taken from the GQ field. If GQ (or the desired
##' field) is missing from both FORMAT or INFO, QUAL will be used.
##'
##' The 'sample.name' argument can be used to select genotypes for specific sample from the VCF. In addition,
##' variants that are homozygous reference in this sample will be filtered. If 'sample.name' is not in
##' the VCF, the first sample will be selected (default). To force the entire VCF to be read no matter the
##' genotypes of samples, use 'sample.name=NULL'.
##'
##' Alleles are split and, for each, column 'ac' reports the allele count. Notable cases incude
##' 'ac=-1' for no/missing calls (e.g. './.'), and 'ac=0' on the first allele to report hom ref,
##' variants. These cases are often filtered later with 'ac>0' to keep only non-ref calls. If
##' the VCF contains no samples or if no sample selection if forced (sample.name=NULL), 'ac' will
##' contain '-1' for all variants in the VCF.
##' @title Read SVs from a VCF file
##' @param vcf.file the path to the VCF file
##' @param keep.ins.seq should it keep the inserted sequence? Default is FALSE.
##' @param keep.ref.seq should it keep the reference allele sequence? Default is FALSE.
##' @param sample.name the name of the sample to use. If "" (default) or sample names not in the VCF,
##' select the first sample. If NULL, don't select particular sample.
##' @param qual.field field to use as quality. Can be in INFO (e.g. default GQ) or
##' FORMAT (e.g. DP). If not found in INFO/FORMAT, QUAL field is used.
##' @param other.field name of another field to extract from the INFO (e.g. AF). Default is NULL
##' @param check.inv should the sequence of MNV be compared to identify inversions. 
##' @param keep.ids keep variant ids? Default is FALSE.
##' @param nocalls if TRUE returns no-calls only (genotype ./.). Default FALSE.
##' @param out.fmt output format. Default is 'gr' for GRanges. Other options: 'df' for
##' data.frame and 'vcf' for the VCF object from the VariantAnnotation package.
##' @param min.sv.size the minimum size of the variant to extract from the VCF. Default is 10 
##' @return depending on 'out.fmt', a GRanges, data.frame, or VCF object with relevant information.
##' @author Jean Monlong
##' @export
##' @examples
##' \dontrun{
##' calls.gr = readSVvcf('calls.vcf')
##' }
readSVvcf <- function(vcf.file, keep.ins.seq=FALSE, keep.ref.seq=FALSE, sample.name='',
                      qual.field=c('GQ', 'QUAL'), other.field=NULL, check.inv=FALSE,
                      keep.ids=FALSE, nocalls=FALSE, out.fmt=c('gr', 'df', 'vcf'),
                      min.sv.size=10){
  ## check file path
  vcf.file = path.expand(vcf.file)
  if(!file.exists(vcf.file)){
    stop(vcf.file, ' not found')
  }

  ## guess if gzipped or not
  con = file(vcf.file)
  use_gz = FALSE
  if(summary(con)$class == 'gzfile'){
    use_gz = TRUE
  }
  close(con)
  
  ## read VCF into a data.frame
  if(is.null(sample.name)){
    ## if NULL, read the entire VCF, encoded with a * in the cpp function
    sample.name = '*'
  }
  if(is.null(other.field)){
    other.field = ''
  }
  svs = read_vcf_cpp(vcf.file, use_gz, sample_name=sample.name,
                     shorten_ref=!keep.ref.seq, shorten_alt=!keep.ins.seq,
                     check_inv=check.inv, gq_field=qual.field[1],
                     keep_nocalls=nocalls, other_field=other.field,
                     min_sv_size=min.sv.size)

  ## Add missing information
  svs$qual = ifelse(svs$qual==-1, NA, svs$qual)
  if(all(is.na(svs$qual))){
    ## needs to be "numeric" for the VCF object output...
    svs$qual = as.numeric(rep(NA, length(svs$qual)))
  }

  ## Eventually convert the additional column extracted
  if(other.field %in% colnames(svs)){
    svs[, other.field] = utils::type.convert(svs[, other.field], as.is=TRUE)
  }

  ## Remove ids if we don't want them
  if(!keep.ids){
    svs$svid = NULL
  }
  
  ## Convert to GRanges or VCF object
  if(out.fmt[1] == 'gr'){
    if(nrow(svs)==0){
      svs = GenomicRanges::GRanges()
    } else {
      svs = GenomicRanges::makeGRangesFromDataFrame(svs, keep.extra.columns=TRUE)
    }
  }
  if(out.fmt[1] == 'vcf'){

    if(nrow(svs)==0){
      vcf.o = VariantAnnotation::VCF()
    } else {
      vcf.o = VariantAnnotation::VCF(GenomicRanges::GRanges(
                                                      svs$seqnames,
                                                      IRanges::IRanges(svs$start, svs$end),
                                                      ),
                                     collapsed=FALSE)

      VariantAnnotation::qual(vcf.o) = svs$qual
      VariantAnnotation::ref(vcf.o) = Biostrings::DNAStringSet(svs$ref)
      VariantAnnotation::alt(vcf.o) = svs$alt

      info.h = S4Vectors::DataFrame(
                            Number=rep('1', 4),
                            Type=c(rep('Integer', 3), 'String'),
                            Description=c(
                              'End coordinate',
                              'SV length (bp)',
                              'Allele count (1:heterozygous, 2:homozygous)',
                              'SV type (DEL, INS, INV)'))
      rownames(info.h) = c('END', 'SVLEN', 'AC', 'SVTYPE')

      info.df = S4Vectors::DataFrame(
                             END=svs$end,
                             SVLEN=svs$size,
                             AC=svs$ac,
                             SVTYPE=svs$type)
      
      ## other field
      if(other.field != '' & other.field %in% colnames(svs)){
        info.n = rownames(info.h)
        info.h = rbind(info.h,
                       x=S4Vectors::DataFrame(Number='1',
                                              Type=ifelse(is.numeric(svs[,other.field]), 'Float', 'String'),
                                              Description=''))
        rownames(info.h) = c(info.n, other.field)
        info.df = cbind(info.df, S4Vectors::DataFrame(OTHER=svs[,other.field]))
        colnames(info.df)[ncol(info.df)] = other.field
      }
      
      VariantAnnotation::info(VariantAnnotation::header(vcf.o)) = info.h
      VariantAnnotation::info(vcf.o) = info.df
      names(vcf.o) = svs$svid

      ## unname geno (otherwise error when using writeVcf)
      names(VariantAnnotation::geno(vcf.o)) = character(0)
    }
    
    svs = vcf.o
  }

  return(svs)  
}
