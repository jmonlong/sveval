% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{read_vcf_multisamps_cpp}
\alias{read_vcf_multisamps_cpp}
\title{Read VCF using CPP reader}
\usage{
read_vcf_multisamps_cpp(
  filename,
  use_gz,
  min_sv_size = 10L,
  shorten_ref = TRUE,
  shorten_alt = TRUE,
  check_inv = FALSE
)
}
\arguments{
\item{filename}{the path to the VCF file (unzipped or gzipped).}

\item{use_gz}{is the VCF file gzipped?}

\item{min_sv_size}{minimum variant size to keep in bp. Variants shorter than this
will be skipped. Default is 10.}

\item{shorten_ref}{should the REF sequence be shortened to the first 10 bp. Default is TRUE}

\item{shorten_alt}{should the ALT sequence be shortened to the first 10 bp. Default is TRUE}

\item{check_inv}{guess if a variant is an inversion by aligning REF with the
reverse complement of ALT. If >80\% similar (and REF and ALT>10bp), variant is classified as INV.}
}
\value{
data.frame with variant and genotype information
}
\description{
For each VCF record the information in the INFO field is used in priority. If
missing, information is guessed from the REF/ALT sequences.
If multiple alleles are defined in ALT, they are split and the allele count extracted
from the GT field.
}
\details{
Alleles are split and, for each, the allele count is computed across samples.
}
\author{
Jean Monlong
}
\keyword{internal}
