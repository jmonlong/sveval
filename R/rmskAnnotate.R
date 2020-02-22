##' Extracts REF/ALT sequence from a VCF, runs RepeatMasker to annotate transposable
##' elements and simple repeats, and annotates the original variants.
##'
##' This is a simple annotation where only the main repeat class is retrieved for
##' each variant (covering most of the sequence).
##'
##' RepeatMasker must be installed.
##' @title Annotate REF/ALT sequence with RepeatMasker
##' @param svs.gr SVs. A GRanges or the path to a VCF file.
##' @param nb.cores the number of cores that RepeatMasker should use. Default:1.
##' @param species the species to use in RepeatMasker. Default: human.
##' @param docker.image docker image with RepeatMasker. Default is NULL.
##' @return an updated GRanges with new columns 'rmsk.name', 'rmsk.classfam' and 'rmsk.cov'.
##' @author Jean Monlong
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @export
##' @examples
##' \dontrun{
##' svs.gr = readSVvcf('calls.vcf', keep.ins.seq=TRUE, keep.ref.seq=TRUE)
##' svs.gr = rmskAnnotate(svs.gr)
##' }
rmskAnnotate <- function(svs.gr, nb.cores=1, species='human', docker.image=NULL){
  if(is.character(svs.gr) & length(svs.gr)==1){
    svs.gr = readSVvcf(svs.gr, keep.ins.seq=TRUE, keep.ref.seq=TRUE)
  }
  svs.gr$id = paste(as.character(GenomicRanges::seqnames(svs.gr)),
                    GenomicRanges::start(svs.gr), svs.gr$type, svs.gr$size, sep='_')
  ## Write sequence FASTA for each variant
  ## seqs = Biostrings::DNAString(ifelse(svs.gr$type %in% c('DEL', 'INV', 'DUP'), svs.gr$REF, svs.gr$ALT))
  seqs = svs.gr$ALT
  seqs[which(svs.gr$type %in% c('DEL', 'INV', 'DUP'))] = svs.gr$REF[which(svs.gr$type %in% c('DEL', 'INV', 'DUP'))]
  names(seqs) = svs.gr$id
  temp.fa = paste0(tempfile(), '.fa')
  Biostrings::writeXStringSet(seqs, temp.fa, format='FASTA')
  ## Run RepeatMasker
  if(!is.null(docker.image)){
    temp.dir = dirname(temp.fa)
    system2('docker', c('run', '-t', '-v', paste0(temp.dir, ':', temp.dir), docker.image, 'RepeatMasker', temp.fa, '--species', species, '-pa', nb.cores))
  } else {
    system2('RepeatMasker', c(temp.fa, '--species', species, '-pa', nb.cores))
  }
  ## Parse RepeatMasker's output
  rmout = utils::read.table(paste0(temp.fa, '.out'), skip=3, as.is=TRUE, fill=TRUE)
  rmout = rmout[,c(5:7, 10,11)]
  colnames(rmout) = c('id', 'start', 'end', 'repeat.name', 'repeat.class.family')
  rmout$rm.w = rmout$end - rmout$start
  rmout = dplyr::group_by(rmout, .data$id) %>% dplyr::arrange(dplyr::desc(.data$rm.w)) %>% dplyr::do(utils::head(.data, 1)) %>%
    as.data.frame(stringsAsFactors=FALSE)
  rownames(rmout) = rmout$id
  ## Annotate SVs
  svs.gr$rmsk.classfam = rmout[svs.gr$id,'repeat.class.family']
  svs.gr$rmsk.name = rmout[svs.gr$id,'repeat.name']
  svs.gr$rmsk.cov = rmout[svs.gr$id,'rm.w'] / GenomicRanges::width(seqs)
  return(svs.gr)
}
