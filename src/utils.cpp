#include <vector>
#include <string>

// split a string
std::vector<std::string> split_str(std::string in_str, std::string delim="\t"){
  std::vector<std::string> line_v;
  auto f_start = 0U;
  auto f_end = in_str.find(delim);
  while (f_end != std::string::npos)
    {
      line_v.push_back(in_str.substr(f_start, f_end - f_start));
      f_start = f_end + delim.length();
      f_end = in_str.find(delim, f_start);
    }
  line_v.push_back(in_str.substr(f_start, in_str.size()));
  return line_v;
}

// reverse complement
std::string rev_comp(std::string seq){
  std::string rc = seq;
  for(unsigned int ii=0; ii<seq.length(); ii++){
    if(seq[ii] == 'A'){
      rc[seq.length()-ii-1] = 'T';
    } else if(seq[ii] == 'T'){
      rc[seq.length()-ii-1] = 'A';
    } else if(seq[ii] == 'G'){
      rc[seq.length()-ii-1] = 'C';
    } else if(seq[ii] == 'C'){
      rc[seq.length()-ii-1] = 'G';
    } 
  }
  return rc;
}
