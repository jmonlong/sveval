##' Wrapper functions that reads arguments and runs sveval.
##' @title Wrapper function to use sveval as a command line tool
##' @param subcmd which subcommand to use. Default is '' which prints the list of available subcommands
##' @param args vector of arguments. By default, reads them from the command line.
##' @return return code (0: success, 1: error)
##' @author Jean Monlong
##' @export
wrapper <- function(subcmd=c('', 'sveval', 'mergetsvs'), args=commandArgs(TRUE)){
  ## print help message about subcommands if no arguments
  printHelp <- function(){
    help.msg = paste0(
      'Subcommands: \n',
      '\tsveval: \n\t\tevaluate SV calls against a truthset\n',
      '\tmergetsvs: \n\t\tmerge TSVs from running sveval on different calls\n\t\t(Warning: files must be named as in the Snakemake pipeline)\n',
      ''
    )
    message(help.msg)
  }

  if(subcmd[1] == '' & length(args) == 0){
    printHelp()
    q(status=1)
  }
  
  ## run the subcommands
  if(subcmd[1] == 'sveval'){
    wrapper_sveval(args)
  } else if (subcmd[1] == 'mergetsvs'){
    wrapper_mergetsvs(args)
  } else {
    stop("Unknown subcommand '", subcmd[1], "'. Must be one of: ", paste(c('sveval', 'mergetsvs'), collapse=' '))
  }
}
