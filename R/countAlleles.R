##' Read a VCF file and count alleles for SVs grouped by SV site. The \code{sv.sites} parameter
##' is a list defining which SVs are assigned to which SV site. The SV ids are those
##' created by \code{readSVvcf.multisamps} which is expected to be run first.
##' @title Count alleles in SV sites across samples
##' @param vcf.file the path to the VCF file
##' @param sv.sites a list defining the SV sites (each element contains a vector with SV ids)
##' @param gq.instead return a matrix with the minimum genotype quality for each site instead of the
##' allele counts. Default is FALSE.
##' @return a matrix with allele counts for each site (rows) and samples (columns)
##' @author Jean Monlong
##' @export
countAlleles <- function(vcf.file, sv.sites, gq.instead=FALSE){
  ## if not named, create names
  if(is.null(names(sv.sites))) {
    names(sv.sites) = paste0('svsite', 1:length(sv.sites))
  }
  if(!is.list(sv.sites)){
    stop('sv.sites is not a list')
  }
  
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
  if(gq.instead){
    return(merge_gq_svsites_cpp(vcf.file, use_gz, sv.sites))
  } else {
    return(merge_ac_svsites_cpp(vcf.file, use_gz, sv.sites))
  }    
}
