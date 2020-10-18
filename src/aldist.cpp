#include <Rcpp.h>
#include "edlib.h" // from https://github.com/Martinsos/edlib
using namespace Rcpp;

//' Using edlib to quickly compute the edit distances between pairs of sequences
//' @title Align sequences and compute edit distance
//' @param seq1 a vector with first sequences in the pair to align
//' @param seq2 a vector with second sequences in the pair to align
//' @return a vector with edit distances
//' @author Jean Monlong
//' @keywords internal
// [[Rcpp::export]]
std::vector<int> aldist(std::vector<std::string> seq1, std::vector<std::string> seq2){
  std::vector<int> res;
  for(unsigned int ii=0; ii<seq1.size(); ii++){
    // compare sequences
    EdlibAlignConfig edconf = edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0);
    EdlibAlignResult ed_o = edlibAlign(seq1[ii].c_str(), seq1[ii].length(),
                                       seq2[ii].c_str(), seq2[ii].length(),
                                       edconf);
    res.push_back(ed_o.editDistance);
  }
  return res;
}
