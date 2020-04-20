// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// read_vcf_cpp
DataFrame read_vcf_cpp(std::string filename, bool use_gz, std::string sample_name, int min_sv_size, bool shorten_ref, bool shorten_alt, std::string gq_field, bool check_inv, bool keep_nocalls, std::string other_field);
RcppExport SEXP _sveval_read_vcf_cpp(SEXP filenameSEXP, SEXP use_gzSEXP, SEXP sample_nameSEXP, SEXP min_sv_sizeSEXP, SEXP shorten_refSEXP, SEXP shorten_altSEXP, SEXP gq_fieldSEXP, SEXP check_invSEXP, SEXP keep_nocallsSEXP, SEXP other_fieldSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< bool >::type use_gz(use_gzSEXP);
    Rcpp::traits::input_parameter< std::string >::type sample_name(sample_nameSEXP);
    Rcpp::traits::input_parameter< int >::type min_sv_size(min_sv_sizeSEXP);
    Rcpp::traits::input_parameter< bool >::type shorten_ref(shorten_refSEXP);
    Rcpp::traits::input_parameter< bool >::type shorten_alt(shorten_altSEXP);
    Rcpp::traits::input_parameter< std::string >::type gq_field(gq_fieldSEXP);
    Rcpp::traits::input_parameter< bool >::type check_inv(check_invSEXP);
    Rcpp::traits::input_parameter< bool >::type keep_nocalls(keep_nocallsSEXP);
    Rcpp::traits::input_parameter< std::string >::type other_field(other_fieldSEXP);
    rcpp_result_gen = Rcpp::wrap(read_vcf_cpp(filename, use_gz, sample_name, min_sv_size, shorten_ref, shorten_alt, gq_field, check_inv, keep_nocalls, other_field));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sveval_read_vcf_cpp", (DL_FUNC) &_sveval_read_vcf_cpp, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_sveval(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
