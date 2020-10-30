#include <Rcpp.h>
#include <fstream>
#include "gzstream.h" // from https://gitlab.com/hrbrmstr/ndjson/-/blob/master/src/gzstream.h
#include "utils.h"
using namespace Rcpp;

// // split a string
// std::vector<std::string> split_str(std::string in_str, std::string delim="\t"){
//   std::vector<std::string> line_v;
//   auto f_start = 0U;
//   auto f_end = in_str.find(delim);
//   while (f_end != std::string::npos)
//     {
//       line_v.push_back(in_str.substr(f_start, f_end - f_start));
//       f_start = f_end + delim.length();
//       f_end = in_str.find(delim, f_start);
//     }
//   line_v.push_back(in_str.substr(f_start, in_str.size()));
//   return line_v;
// }

//' Count allele at the SV site level. This function will use variant ids
//' created by read_vcf_multisamps_cpp
//' @title Count alleles in SV sites across samples
//' @param filename the path to the VCF file (unzipped or gzipped).
//' @param use_gz is the VCF file gzipped?
//' @param sv_sites a list defining the variant ID for each sv site (element). 
//' @return matrix with allele counts for each sv site (rows) and sample (columns)
//' @author Jean Monlong
//' @keywords internal
// [[Rcpp::export]]
IntegerMatrix merge_ac_svsites_cpp(std::string filename, bool use_gz, List sv_sites){
  // make hash to associate a variant id to a sv site
  std::map<std::string, int> sv_to_site_i;
  std::vector<std::string> site_ids = sv_sites.attr("names");
  for(int ii=0; ii < sv_sites.size(); ii++){
    std::string site_id = site_ids[ii];
    std::vector<std::string> svs = as<std::vector<std::string>>(sv_sites[ii]);
    for(unsigned int jj=0; jj < svs.size(); jj++){
      sv_to_site_i[svs[jj]] = ii;
    }    
  }

  // set up file reader
  std::ifstream inh_file;
  igzstream inh_file_gz;
  std::string lineh;
    
  // read beginning of VCF once to list/count samples
  std::vector<std::string> samps;
  if(use_gz){
    inh_file_gz.open(filename.c_str());
  } else {
    inh_file.open(filename);
  }
  bool getmore = false;
  if(use_gz){
    getmore = bool(std::getline(inh_file_gz, lineh));
  } else {
    getmore = bool(std::getline(inh_file, lineh));
  }
  while (getmore) {
    Rcpp::checkUserInterrupt();
    // skip if header lines
    if(lineh[0] == '#'){
      if(lineh[1] != '#'){
        // header line with column names including sample names
        std::vector<std::string> header_v = split_str(lineh);
        if((header_v.size() >= 10)){
          for(unsigned int ss=9; ss<header_v.size(); ss++){
            samps.push_back(header_v[ss]);
          }
          getmore = false;
          continue;
        }
      }
    }
    if(use_gz){
      getmore = bool(std::getline(inh_file_gz, lineh));
    } else {
      getmore = bool(std::getline(inh_file, lineh));
    }
  }
  inh_file.close();
  inh_file_gz.close();

  // init output matrix
  IntegerMatrix ac_mat(site_ids.size(), samps.size());
  colnames(ac_mat) = wrap(samps);
  rownames(ac_mat) = wrap(site_ids);

  // read the VCF entirely this time to update the matrix
  std::ifstream in_file;
  igzstream in_file_gz;
  std::string line;
  if(use_gz){
    in_file_gz.open(filename.c_str());
  } else {
    in_file.open(filename);
  }
  if(use_gz){
    getmore = bool(std::getline(in_file_gz, line));
  } else {
    getmore = bool(std::getline(in_file, line));
  }
  int line_id = 0;
  while (getmore) {
    Rcpp::checkUserInterrupt();
    // skip if header lines
    if(line[0] == '#'){
      if(use_gz){
        getmore = bool(std::getline(in_file_gz, line));
      } else {
        getmore = bool(std::getline(in_file, line));
      }
      continue;
    }
    
    // otherwise tab-split line
    std::vector<std::string> line_v = split_str(line);
    // std::cout << "line size: " << line_v.size() << std::endl;

    // parse INFO
    std::vector<std::string> infos_v = split_str(line_v[7], ";");
    std::map<std::string,std::string> infos;
    for (unsigned int ii=0; ii<infos_v.size(); ii++){
      std::vector<std::string> info_pair = split_str(infos_v[ii], "=");
      if(info_pair.size() == 2){
        infos[info_pair[0]] = info_pair[1];
      } else {
        infos[infos_v[ii]] = "";
      }
    }
        
    //
    // parse the alleles: split multi-ALT and update allele count matrix 
    //
    // check if any allele in the set of interest. if not, might as well skip the variant
    int nb_alts = split_str(line_v[4], ",").size();
    bool in_svsites = false;
    for(int al_id=0; al_id<nb_alts; al_id++){
      std::string svid = "sv_" + std::to_string(line_id) + "_" + std::to_string(al_id);
      if(sv_to_site_i.find(svid) != sv_to_site_i.end()){
        in_svsites = true;
      }
    }
    // if at least one allele is in the set of interest
    if(in_svsites){
      // parse genotype for each sample
      std::vector<std::string> format_fields_v = split_str(line_v[8], ":");
      for(unsigned int sample_col=9; sample_col<samps.size()+9; sample_col++){
        // parse FORMAT/GT fields
        std::vector<std::string> gt_fields_v = split_str(line_v[sample_col], ":");
        std::map<std::string,std::string> gt_fields;
        for (unsigned int ii=0; ii<gt_fields_v.size(); ii++){
          gt_fields[format_fields_v[ii]] = gt_fields_v[ii];
        }
        // split GT
        if(gt_fields.find("GT") == gt_fields.end()){
          Rcout << "GT is missing from FORMAT and genotype field, exiting." << std::endl;
        } else {
          std::string gt_value = gt_fields["GT"];
          std::vector<std::string> gt_s;
          if(gt_value.find("/") != std::string::npos){
            gt_s = split_str(gt_value, "/");
          } else if (gt_value.find("|") != std::string::npos){
            gt_s = split_str(gt_value, "|");
          } else {
            gt_s.push_back(gt_value);
          }
          // update matrix for every non-ref allele found
          for(unsigned int ii=0; ii<gt_s.size(); ii++){
            if(gt_s[ii] != "." && gt_s[ii] != "0"){
              int al_id = atoi(gt_s[ii].c_str()) - 1;
              std::string svid = "sv_" + std::to_string(line_id) + "_" + std::to_string(al_id);
              // std::string samp = samps[sample_col-9];
              if(sv_to_site_i.find(svid) != sv_to_site_i.end()){
                ac_mat(sv_to_site_i[svid], sample_col-9) = ac_mat(sv_to_site_i[svid], sample_col-9) + 1;
              }
            }
          }
        }
      }
    }
    Rcpp::checkUserInterrupt(); 
    line_id++;
    // get next line
    if(use_gz){
      getmore = bool(std::getline(in_file_gz, line));
    } else {
      getmore = bool(std::getline(in_file, line));
    }
  }
  in_file.close();
  in_file_gz.close();
  
  // return matrix
  return(ac_mat);
}
