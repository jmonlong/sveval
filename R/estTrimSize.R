##' @title Estimate the amount of sequence to right-trim
##' @param seqs1 a set of sequences, e.g. the reference sequences.
##' @param seqs2 a set of sequences, e.g. the alt sequences.
##' @return the number of bases to right trim.
##' @author Jean Monlong
##' @keywords internal
estTrimSize <- function(seqs1, seqs2){
  ## Get longest common suffix
  res = sapply(1:length(seqs1), function(ii) Biostrings::lcsuffix(seqs1[[ii]], seqs2[[ii]]))
  ## If same size as one of the sequence, make it one base shorter
  res = ifelse(res==Biostrings::nchar(seqs1), res-1, res)
  res = ifelse(res==Biostrings::nchar(seqs2), res-1, res)
  return(res)
}
