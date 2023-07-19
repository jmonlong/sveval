##' Wrapper functions that reads arguments and runs sveval.
##' @title Internal wrapper function to run sveval as a command line tool
##' @param args arguments
##' @return return code (0: success, 1: error)
##' @author Jean Monlong
##' @keywords internal
wrapper_sveval <- function(args){
  printHelp <- function(args, do.print=TRUE, desc=NULL){
    help.msg = lapply(args, function(aa){
      paste0('\t', paste(aa$args, collapse=', '), '\t', aa$desc)
    })
    help.msg = paste(unlist(help.msg), collapse='\n')
    if(do.print){
      if(!is.null(desc)){
        message(desc)
      }
      message(help.msg)
    }
    return(help.msg)
  }
  
  args = commandArgs(TRUE)

  ## define arguments (no "value" means required)
  cmd.args = list()
  cmd.args$calls = list(desc="input VCF with the SV calls", args=c("-c", "-calls"))
  cmd.args$truth = list(desc="input VCF with the truth SVs", args=c("-t", "-truth"))
  cmd.args$sample = list(desc="sample name", args=c("-s", "-sample"))
  cmd.args$output = list(desc="output prefix for Rdata and TSV files", args=c("-o", "-output"))
  cmd.args$eval = list(value='call', desc="Optional. evaluation type: 'call' for absence/presence, 'geno' if genotypes have to match. Default: call", args=c("-e", "-eval"))
  cmd.args$region = list(value=NA, desc="Optional. regions of interest. NA means whole genome. Default: NA", args=c("-r", "-region"))
  cmd.args$simprep = list(value=NA, desc="Optional. simple repeat track. NA means disabled. Default: NA", args=c("-p", "-simprep"))
  cmd.args$inversion = list(value=FALSE, desc="Optional. look for inversions? Default: FALSE", args=c("-i", "-inversion"))
  cmd.args$minol = list(value=0.5, desc="Optional. minimum overlap to match variants. Default:0.5", args=c("-m", "-minol"))

  ## print help message if no arguments or help arguments
  if(length(args)==0 || args[1] %in% c('-h', '-help')){
    printHelp(cmd.args, desc='Evaluate SV calls against a truthset')
    q(status=1)
  }

  ## otherwise parse arguments
  arg.ii = 1
  while(arg.ii < length(args)){
    for(argn in names(cmd.args)){
      if(args[arg.ii] %in% cmd.args[[argn]][["args"]]){
        cmd.args[[argn]][["value"]] = args[arg.ii+1]
      }
    }
    arg.ii = arg.ii + 2
  }

  ## check that all arguments have a value
  missing.args = which(unlist(lapply(cmd.args, function(x) all(names(x)!='value'))))
  if(length(missing.args)>0){
    stop("Missing values for arguments:\n", printHelp(cmd.args[missing.args], do.print=FALSE))
  }
  
  ## Regions of interest or whole-genome?
  bed = NULL
  if(cmd.args$region$value != 'NA'){
    bed = cmd.args$region$value
  }

  ## simple repeat track
  sr = NULL
  if(cmd.args$simprep$value != 'NA'){
    sr = utils::read.table(cmd.args$simprep$value, as.is=TRUE, sep='\t')
    sr = GenomicRanges::GRanges(sr[,1], IRanges::IRanges(sr[,2], sr[,3]))
  }

  ## evaluation
  eval.o = svevalOl(cmd.args$calls$value, cmd.args$truth$value,
                    sample.name=cmd.args$sample$value,
                    check.inv=as.logical(cmd.args$inversion$value),
                    max.ins.dist=100,
                    bed.regions=bed,
                    min.ol=as.numeric(cmd.args$minol$value), 
                    geno.eval=as.logical(cmd.args$eval$value=='geno'),
                    stitch.hets=as.logical(cmd.args$eval$value=='geno'),
                    merge.hets=as.logical(cmd.args$eval$value=='geno'),
                    simprep=sr,
                    min.size=50)

  ## save RData object
  save(eval.o, file=paste0(cmd.args$output$value, '-sveval.RData'))

  ## PR curve
  utils::write.table(eval.o$curve, file=paste0(cmd.args$output$value, '-prcurve.tsv'),
              sep='\t', quote=FALSE, row.names=FALSE)
  
  ## per-size estimates
  df = plot_persize(eval.o, plot=FALSE,
                    size.breaks=c(50,100,200,300,400,600,800,
                                  1000,2500,5000,Inf))
  utils::write.table(df, file=paste0(cmd.args$output$value, '-persize.tsv'),
              sep='\t', quote=FALSE, row.names=FALSE)

  q(status=0)

}
