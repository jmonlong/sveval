% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepareOl.R
\name{prepareOl}
\alias{prepareOl}
\title{Prepare SV overlaps before annotation}
\usage{
prepareOl(
  query,
  subject,
  min.rol = 0.1,
  max.ins.dist = 1,
  range.seq.comp = FALSE,
  ins.seq.comp = FALSE,
  nb.cores = 1,
  by.gt = FALSE
)
}
\arguments{
\item{query}{a query GRanges object}

\item{subject}{a subject GRanges object}

\item{min.rol}{minimum reciprocal overlap for deletions and other "ranges" SVs. Default is 0.1}

\item{max.ins.dist}{maximum distance for insertions to be clustered.}

\item{range.seq.comp}{compare sequence instead of overlapping deletions/inversion/etc. Default is FALSE.}

\item{ins.seq.comp}{compare sequence instead of insertion sizes. Default is FALSE.}

\item{nb.cores}{number of processors to use. Default is 1.}

\item{by.gt}{should the variants be split by genotype? Default is FALSE, i.e. all variants with
an alternate allele (ac>0) is considered 'called'.}
}
\value{
a GRanges with information about pairs of SVs in query and subject that overlap
\item{GRange}{intersected ranges (for "ranges" SVs)}
\item{queryHits}{the id of the input query}
\item{subjectHits}{the id of the input subject}
\item{querSize}{the size of the input query}
\item{subjectSize}{the size of the input subject}
\item{interSize}{the size of the intersection (e.g. range, ins size, ins seq alignment)}
\item{type}{the SV type of the pair}
}
\description{
Prepare SV overlaps before annotation
}
\author{
Jean Monlong
}
\keyword{internal}
