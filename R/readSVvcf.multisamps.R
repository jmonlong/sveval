##' Read a VCF file that contains SVs for multiple samples and create a GRanges with population estimates.
##' 
##' Alleles are split and, for each, column 'ac' reports the total allele count. The number of "ref" genotypes,
##' i.e. homozygous refs e.g. '0/0', and the number of called samples (not missing './.') are also counted.
##' The 'af' column contains an estimate of the allele frequency.
##' @title Read SVs from a VCF file
##' @param vcf.file the path to the VCF file
##' @param keep.ins.seq should it keep the inserted sequence? Default is FALSE.
##' @param keep.ref.seq should it keep the reference allele sequence? Default is FALSE.
##' @param check.inv should the sequence of MNV be compared to identify inversions. 
##' @param keep.ids keep variant ids? Default is FALSE.
##' @param out.fmt output format. Default is 'gr' for GRanges. Other options: 'df' for
##' data.frame.
##' @param min.sv.size the minimum size of the variant to extract from the VCF. Default is 10 
##' @return depending on 'out.fmt', a GRanges or a data.frame with relevant information.
##' @author Jean Monlong
##' @export
##' @examples
##' \dontrun{
##' svs.gr = readSVvcf.multisamps('svs.vcf')
##' }
readSVvcf.multisamps <- function(vcf.file, keep.ins.seq=FALSE, keep.ref.seq=FALSE,
                                 check.inv=FALSE, keep.ids=FALSE,
                                 out.fmt=c('gr', 'df'), min.sv.size=10){
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
  svs = read_vcf_multisamps_cpp(vcf.file, use_gz, min_sv_size=min.sv.size,
                                shorten_ref=!keep.ref.seq, shorten_alt=!keep.ins.seq,
                                check_inv=check.inv)

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

  return(svs)  
}
