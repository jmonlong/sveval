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
//' Alleles are split and, for each, column 'ac' reports the allele count. Notable cases incude
//' 'ac=-1' for no/missing calls (e.g. './.'), and 'ac=0' on the first allele to report hom ref,
//' variants. These cases are often filtered later with 'ac>0' to keep only non-ref calls. If
//' the VCF contains no samples or if no sample selection if forced (sample_name='*'), 'ac' will
//' contain '-1' for all variants in the VCF.
//' @title Read VCF using CPP reader
//' @param filename the path to the VCF file (unzipped or gzipped).
//' @param use_gz is the VCF file gzipped?
//' @param sample_name which sample to process. If not found, uses first sample in VCF file.
//' If "*", force no sample selection
//' @param min_sv_size minimum variant size to keep in bp. Variants shorter than this
//' will be skipped. Default is 10. 
//' @param shorten_ref should the REF sequence be shortened to the first 10 bp. Default is TRUE
//' @param shorten_alt should the ALT sequence be shortened to the first 10 bp. Default is TRUE
//' @param gq_field which field from FORMAT should be used as genotype quality. Default is "GQ".
//' If not found, QUAL will be used
//' @param check_inv guess if a variant is an inversion by aligning REF with the
//' reverse complement of ALT. If >80\% similar (and REF and ALT>10bp), variant is classified as INV.
//' @param keep_nocalls should we keep variants/alleles with missing genotypes (e.g. "./.").
//' Default is FALSE
//' @param other_fields name of another field from INFO to extract.
//' @return data.frame with variant and genotype information
//' @author Jean Monlong
//' @keywords internal
// [[Rcpp::export]]
DataFrame read_vcf_cpp(std::string filename, bool use_gz, std::string sample_name="", int min_sv_size=10,
                       bool shorten_ref=true, bool shorten_alt=true, std::string gq_field="GQ",
                       bool check_inv=false, bool keep_nocalls=false, CharacterVector other_fields=CharacterVector::create()){
  // info to extract. will be the columns of the output dataframe
  std::vector<std::string> seqnames;                      // chromosome name
  std::vector<int> starts;                                // start position
  std::vector<int> ends;                                  // end position
  std::vector<std::string> svids;                         // variant ID
  std::vector<std::string> filters;                       // FILTER values
  std::vector<std::string> refs;                          // reference allele sequence
  std::vector<std::string> alts;                          // alternate allele sequence
  std::vector<std::string> svtypes;                       // SV type (INS, DEL, INV, ...)
  std::vector<int> sizes;                                 // SV size (absolute value)
  std::vector<int> acs;                                   // allele counts (for a diploid genome: 1=heterozygous, 2=homozygous)
  std::vector<double> quals;                              // genotype quality

  // init object for potential other fields to save
  bool found_other=false; // was the additional found, aka should it be added to the output
  std::map<std::string,std::vector<std::string>> others;  // other columns to extract
  for(CharacterVector::iterator it=other_fields.begin(); it != other_fields.end(); it++){
    String s = *it; // then use s.get_cstring()
    std::vector<std::string> newcol;
    others[s] = newcol;
  }
  
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
  int sample_col=-1;

  while (getmore) {
    Rcpp::checkUserInterrupt();
    // skip if header lines
    if(line[0] == '#'){
      if(line[1] != '#'){
        // header line with column names including sample names
        std::vector<std::string> header_v = split_str(line);
        if((header_v.size() >= 10) & (sample_name != "*")){ // samples in VCF
          sample_col = 9; // default is first sample column
          for(unsigned int ii=9; ii<header_v.size(); ii++){
            if(header_v[ii] == sample_name){
              sample_col = ii;
            }
          }
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
    if((size == -1) & (end_rec != -1) & (svtype != "") & (svtype != "INS")){ // size from start and end
      size = end_rec - start_rec;
    }
    
    // init quality with QUAL
    double qual = -1;
    if(line_v[5] != "."){
      qual = atof(line_v[5].c_str());
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
    bool any_nonref = false; // any non-ref allele. if not, include hom ref in output
    if(sample_col == -1){
      // keep all alleles
      for(unsigned int ii=0; ii < alt_seqs.size(); ii++){
        als_count[ii] = -1;
      }
    } else {
      // parse field
      std::vector<std::string> gt_fields_v = split_str(line_v[sample_col], ":");
      std::vector<std::string> format_fields_v = split_str(line_v[8], ":");
      std::map<std::string,std::string> gt_fields;
      for (unsigned int ii=0; ii<gt_fields_v.size(); ii++){
        gt_fields[format_fields_v[ii]] = gt_fields_v[ii];
      }

      // split GT
      if(gt_fields.find("GT") == gt_fields.end()){
        Rcout << "GT is missing from FORMAT and genotype field, exiting." << std::endl;
        als_count[0] = -1;
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
          if(gt_s[ii] == "."){
            // if missing genotype, save the allele but don't add to allele count
            if(als_count.count(0) == 0){
              als_count[0] = -1;
            }
          } else {
            // otherwise increment the allele count of the appropriate allele
            int al_id = atoi(gt_s[ii].c_str()) - 1;
            if(al_id > -1){
              any_nonref = true;
            }
            if(als_count.count(al_id) == 0){
              als_count[al_id] = 1;
            } else {
              als_count[al_id]++;
            }
          }
        }
      }
      
      // find genotype quality from field gq_field
      if(gt_fields.find(gq_field) != gt_fields.end()){
        if(gt_fields[gq_field] == "."){
          qual = -1;
        } else {
          qual = atof(gt_fields[gq_field].c_str());
        }
      } else if (infos.count(gq_field) > 0){
        if(infos[gq_field] == "."){
          qual = -1;
        } else {
          qual = atof(infos[gq_field].c_str());
        }
      }
    }
    
    // for each allele to process, update info and output columns
    for(std::map<int,int>::iterator iter=als_count.begin(); iter!=als_count.end(); iter++){
      int al_id = iter->first;
      int al_count = iter->second;
      if(al_id == -1){
        if(any_nonref){
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
      // skip if too small but keep NBD and translocations
      if((size_al < min_sv_size) & (svtype_al != "BND") & (svtype_al != "TRA")){
        continue;
      }
      // skip if not called, e.g. './.'
      // (except if we want to keep no-calls or we are not selecting for a specific sample)
      if((al_count == -1) & !keep_nocalls & (sample_col > 0)){
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

      // extract additional info field if needed
      for(CharacterVector::iterator it=other_fields.begin(); it != other_fields.end(); it++){
	String other_field = *it;
	if(infos.count(other_field) > 0){
	  others[other_field].push_back(infos[other_field]);
	  found_other = true;
	} else {
	  others[other_field].push_back("");
	}
      }

      // update vectors
      seqnames.push_back(line_v[0]);
      starts.push_back(start_rec);
      ends.push_back(end_rec_al);
      svids.push_back(line_v[2]);
      filters.push_back(line_v[6]);
      acs.push_back(al_count);
      sizes.push_back(size_al);
      svtypes.push_back(svtype_al);
      quals.push_back(qual);
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
			  _["filter"] = filters,
			  _["type"] = svtypes,
			  _["size"] = sizes,
			  _["ac"] = acs,
			  _["ref"] = refs,
			  _["alt"] = alts,
			  _["qual"] = quals,
			  _["stringsAsFactors"] = false);

  if(found_other){
    for(CharacterVector::iterator it=other_fields.begin(); it != other_fields.end(); it++){
      String s = *it; // then use s.get_cstring()
      res.push_back(others[s], s);
    }
  }

  return res;
}
