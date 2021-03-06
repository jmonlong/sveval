#include <Rcpp.h>
#include <fstream>
#include "gzstream.h" // from https://gitlab.com/hrbrmstr/ndjson/-/blob/master/src/gzstream.h
#include "edlib.h" // from https://github.com/Martinsos/edlib
#include "utils.h"
using namespace Rcpp;

//' For each VCF record the information in the INFO field is used in priority. If
//' missing, information is guessed from the REF/ALT sequences.
//' If multiple alleles are defined in ALT, they are split and the allele count extracted
//' from the GT field.
//'
//' Alleles are split and, for each, the allele count is computed across samples. 
//' @title Read VCF using CPP reader
//' @param filename the path to the VCF file (unzipped or gzipped).
//' @param use_gz is the VCF file gzipped?
//' @param min_sv_size minimum variant size to keep in bp. Variants shorter than this
//' will be skipped. Default is 10. 
//' @param shorten_ref should the REF sequence be shortened to the first 10 bp. Default is TRUE
//' @param shorten_alt should the ALT sequence be shortened to the first 10 bp. Default is TRUE
//' @param check_inv guess if a variant is an inversion by aligning REF with the
//' reverse complement of ALT. If >80\% similar (and REF and ALT>10bp), variant is classified as INV.
//' @return data.frame with variant and genotype information
//' @author Jean Monlong
//' @keywords internal
// [[Rcpp::export]]
DataFrame read_vcf_multisamps_cpp(std::string filename, bool use_gz, int min_sv_size=10,
                                  bool shorten_ref=true, bool shorten_alt=true, bool check_inv=false){
  // info to extract. will be the columns of the output dataframe
  std::vector<std::string> seqnames;           // chromosome name
  std::vector<int> starts;                     // start position
  std::vector<int> ends;                       // end position
  std::vector<std::string> svids;              // variant ID
  std::vector<std::string> refs;               // reference allele sequence
  std::vector<std::string> alts;               // alternate allele sequence
  std::vector<std::string> svtypes;            // SV type (INS, DEL, INV, ...)
  std::vector<int> sizes;                      // SV size (absolute value)
  std::vector<int> ac_tots;                    // total allele counts
  std::vector<double> afs;                        // allele frequency
  std::vector<int> nrefs;                      // nb hom refs samples
  std::vector<int> ncalls;                       // nb of sample with any call (!=./.)
  // read file line by line
  std::ifstream in_file;
  igzstream in_file_gz;
  if(use_gz){
    in_file_gz.open(filename.c_str());
  } else {
    in_file.open(filename);
  }
  std::string line;
  bool getmore = false;
  if(use_gz){
    getmore = bool(std::getline(in_file_gz, line));
  } else {
    getmore = bool(std::getline(in_file, line));
  }
  int sample_col_s = -1;
  int sample_col_e = -1;
  int line_id = 0; // for variant id (header line not counted)

  while (getmore) {
    Rcpp::checkUserInterrupt();
    // skip if header lines
    if(line[0] == '#'){
      if(line[1] != '#'){
        // header line with column names including sample names
        std::vector<std::string> header_v = split_str(line);
        if((header_v.size() >= 10)){
          sample_col_s = 9;
          sample_col_e = header_v.size();
        }
      }
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

    // define START coordinate
    int start_rec = atoi(line_v[1].c_str());

    // init size, end, type with INFO if defined
    std::string svtype = "";
    if (infos.count("SVTYPE")>0){ // otherwise use START + SVLEN
      svtype = infos["SVTYPE"];
    }
    int end_rec = -1;
    int size = -1;
    if(infos.count("END")>0){ // use END if there
      end_rec = atoi(infos["END"].c_str());
    }
    if (infos.count("SVLEN")>0){ // SVLEN first
      size = std::abs(atoi(infos["SVLEN"].c_str()));
    } else if (infos.count("INSLEN")>0){ // then INSLEN
      size = atoi(infos["INSLEN"].c_str());
    }
    if ((end_rec == -1) & (size != -1)){ // END from position + size
      if((svtype != "INS") & (svtype != "")){
        end_rec = start_rec + size;
      }
    }
    if((size == -1) & (end_rec != -1)){ // size from start and end
      size = end_rec - start_rec;
    }
        
    //
    // parse the alleles: split multi-ALT, extract quality and
    //    eventually compare REF/ALT sequences to define end/size/type
    
    // split ALT sequences
    std::vector<std::string> alt_seqs = split_str(line_v[4], ",");
    int ref_l = line_v[3].length();

    // rare situation where we don't want to guess the size/end from REF/ALT sequence
    // but some information is missing from the INFO field
    if ((svtype == "INV") & (size == -1)){
      size = ref_l;
    } 

    // decide which alleles to process
    std::map<int,int> als_count;
    int refs_count = 0;
    int calls_count = 0;
    std::vector<std::string> format_fields_v = split_str(line_v[8], ":");
    for(int sample_col=sample_col_s; sample_col<sample_col_e; sample_col++){
      // count allele in one sample
      bool samp_called = false;
      bool hom_ref_samp = true;
      
      // parse field
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
        // count how many times each allele is present.
        // for diploid: 1=het, 2=hom (e.g. 1/1)
        for(unsigned int ii=0; ii<gt_s.size(); ii++){
          if(gt_s[ii] != "."){
            samp_called = true;
            // increment the allele count of the appropriate allele
            int al_id = atoi(gt_s[ii].c_str()) - 1;
            if(al_id >= 0){
              hom_ref_samp = false;
            }
            if(als_count.count(al_id) == 0){
              als_count[al_id] = 1;
            } else {
              als_count[al_id]++;
            }
          }
        }
      }
      if(samp_called){
        calls_count++;
      }
      if(hom_ref_samp){
        refs_count++;
      }
    }
    
    // for each allele to process, update info and output columns
    for(std::map<int,int>::iterator iter=als_count.begin(); iter!=als_count.end(); iter++){
      int al_id = iter->first;
      int al_count = iter->second;
      std::string svid = "sv_" + std::to_string(line_id) + "_" + std::to_string(al_id);
      if(al_id == -1){
        if(als_count.size() > 1){
          // skip ref allele because we found at least one alt allele
          continue;
        } else {
          // output hom ref, i.e. first allele with ac=0
          al_id = 0;
          al_count = 0;
        }
      }
      int alt_l = alt_seqs[al_id].length();
      // update size if necessary
      int size_al = size;
      if (size_al == -1){
        if(alt_l > ref_l){
          size_al = alt_l - ref_l;
        } else {
          size_al = ref_l - alt_l;
        }
      }
      // update type if necessary
      std::string svtype_al = svtype;
      if(svtype_al == ""){
        if(alt_l > ref_l){
          svtype_al = "INS";
        } else {
          svtype_al = "DEL";
        }
        // compare REF vs ALT to guess invesions
        if(check_inv & (alt_l > 10) & (ref_l > 10)){
          std::string alt_rc = rev_comp(alt_seqs[al_id]);
          // compare sequences
          EdlibAlignResult ed_o = edlibAlign(line_v[3].c_str(), ref_l,
                                             alt_rc.c_str(), alt_l,
                                             edlibDefaultAlignConfig());
          int ed = ed_o.editDistance;
          if((ed < 0.2 * alt_l) & (ed < 0.2 * ref_l)){
            svtype_al = "INV";
            size_al = ref_l;
          }
        }
      }
      // skip if too small
      if(size_al < min_sv_size){
        continue;
      }
      // update end if necessary
      int end_rec_al = end_rec;
      if(end_rec_al == -1){
        if(svtype_al == "INS"){
          // simple definition to avoid being too conservative on complex insertions
          end_rec_al = start_rec + ref_l - 1;
        } else {
          end_rec_al = start_rec + size_al;
        }
      }
      
      // update vectors
      seqnames.push_back(line_v[0]);
      starts.push_back(start_rec);
      ends.push_back(end_rec_al);
      svids.push_back(svid);
      ac_tots.push_back(al_count);
      nrefs.push_back(refs_count);
      ncalls.push_back(calls_count);
      // allele frequency as allele count divided by 2 times the samples with calls
      double af = double(al_count) / (2 * calls_count);
      afs.push_back(af);
      sizes.push_back(size_al);
      svtypes.push_back(svtype_al);
      // shorten REF is necessary
      std::string ref_seq = line_v[3];
      if(shorten_ref & (ref_seq.length() > 10)){
        ref_seq = ref_seq.substr(0, 10).insert(10, "...");
      }
      refs.push_back(ref_seq);
      // shorten REF is necessary
      std::string alt_seq = alt_seqs[al_id];
      if(shorten_alt & (alt_seq.length() > 10)){
        alt_seq = alt_seq.substr(0, 10).insert(10, "...");
      }
      alts.push_back(alt_seq);
      
      Rcpp::checkUserInterrupt(); 
    }

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
  
  DataFrame res;
  res = DataFrame::create(
                          _["seqnames"] = seqnames,
                          _["start"] = starts,
                          _["end"] = ends,
                          _["svid"] = svids,
                          _["type"] = svtypes,
                          _["size"] = sizes,
                          _["af"] = afs,
                          _["ac"] = ac_tots,
                          _["nrefs"] = nrefs,
                          _["ncalls"] = ncalls,
                          _["ref"] = refs,
                          _["alt"] = alts,
                          _["stringsAsFactors"] = false);
  return res;
}
