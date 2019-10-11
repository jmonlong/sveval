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
rmskAnnotate <- function(svs.gr, nb.cores=1, species='human'){
  if(is.character(svs.gr) & length(svs.gr)==1){
    svs.gr = readSVvcf(svs.gr, keep.ins.seq=TRUE, keep.ref.seq=TRUE)
  }
  svs.gr$id = paste(as.character(GenomicRanges::seqnames(svs.gr)),
                    GenomicRanges::start(svs.gr), svs.gr$type, svs.gr$size, collapse='_')
  ## Write sequence FASTA for each variant
  seqs = Biostrings::DNAString(ifelse(svs.gr$type %in% c('DEL', 'INV', 'DUP'), svs.gr$REF, svs.gr$ALT))
  names(seqs) = svs.gr$id
  temp.fa = paste0(tempfile(), '.fa')
  Biostrings::writeXStringSet(seqs, temp.fa, format='FASTA')
  ## Run RepeatMasker
  system2('RepeatMasker', c(temp.fa, '--species', species, '-pa', nb.cores))
  ## Parse RepeatMasker's output
  rmout = utils::read.table(paste0(temp.fa, '.out'), skip=3, as.is=TRUE)
  rmout = rmout[,c(5:7, 10,11)]
  colnames(rmout) = c('id', 'start', 'end', 'repeat.name', 'repeat.class.family')
  rmout$rm.w = rmout$end - rmout$start
  rmout = dplyr::group_by(rmout, .data$id) %>% dplyr::arrange(dplyr::desc(.data$rm.w)) %>% dplyr::do(utils::head(.data, 1)) %>%
    as.data.frame(stringsAsFactors=FALSE)
  rownames(rmout) = rmout$id
  ## Annotate SVs
  svs.gr$rmsk.classfam = rmout$repeat.class.family[svs.gr$id]
  svs.gr$rmsk.name = rmout$repeat.name[svs.gr$id]
  svs.gr$rmsk.cov = rmout$rm.w[svs.gr$id] / GenomicRanges::width(seqs)
  return(svs.gr)
}
