% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapper.R
\name{wrapper}
\alias{wrapper}
\title{Wrapper function to use sveval as a command line tool}
\usage{
wrapper(subcmd = c("", "sveval", "mergetsvs"), args = commandArgs(TRUE))
}
\arguments{
\item{subcmd}{which subcommand to use. Default is '' which prints the list of available subcommands}

\item{args}{vector of arguments. By default, reads them from the command line.}
}
\value{
return code (0: success, 1: error)
}
\description{
Wrapper functions that reads arguments and runs sveval.
}
\author{
Jean Monlong
}
