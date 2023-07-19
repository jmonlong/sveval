##' Wrapper functions that merges TSV files made by sveval and make plots
##' @title Internal wrapper function to merge TSVs as a command line tool
##' @param args arguments
##' @return return code (0: success, 1: error)
##' @author Jean Monlong
##' @import ggplot2
##' @importFrom magrittr %>%
##' @importFrom rlang .data
##' @keywords internal
wrapper_mergetsvs <- function(args){
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
  cmd.args$out = list(value="merged-sveval-results",
                      desc="Optional. Prefix for the PDF and TSV outputs. Default: merged-sveval-results",
                      args=c("-o", "-output"))
  cmd.args$tsvs = list(desc="Text file listing all the TSV files for both the PR curve (*-prcurve.tsv) and per-size (*-persize.tsv) evaluation results.",
                      args=c("-t", "-tsvs"))

  ## print help message if no arguments or help arguments
  if(length(args)==0 || args[1] %in% c('-h', '-help')){
    printHelp(cmd.args, desc='Merge TSVs from running sveval on different calls (Warning: files must be named as in the Snakemake pipeline)')
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

  ## read file names for all TSVs
  all.tsvs = scan(cmd.args$tsvs$value, '', quiet=TRUE)

  ## Merge persize files
  tsvs = grep('persize', all.tsvs, value=TRUE)
  df = lapply(tsvs, function(ff){
    df = utils::read.table(ff, as.is=TRUE, header=TRUE)
    ## extract information from the filename (consistent with the Snakemake pipeline)
    df$exp = gsub('(.*)-.*-.*-.*-.*-persize.tsv', '\\1', basename(ff))
    df$method = gsub('.*-(.*)-.*-.*-.*-persize.tsv', '\\1', basename(ff))
    df$sample = gsub('.*-.*-(.*)-.*-.*-persize.tsv', '\\1', basename(ff))
    df$region = gsub('.*-.*-.*-(.*)-.*-persize.tsv', '\\1', basename(ff))
    df$eval = gsub('.*-.*-.*-.*-(.*)-persize.tsv', '\\1', basename(ff))
    df$min.cov = .5
    df
  })
  df = do.call(rbind, df)
  utils::write.table(df, file=paste0(cmd.args$out$value, '-persize.tsv'), quote=FALSE, sep='\t', row.names=FALSE)

  ## Merge PR curve files
  tsvs = grep('prcurve', all.tsvs, value=TRUE)
  df = lapply(tsvs, function(ff){
    df = utils::read.table(ff, as.is=TRUE, header=TRUE)
    ## extract information from the filename
    df$exp = gsub('(.*)-.*-.*-.*-.*-prcurve.tsv', '\\1', basename(ff))
    df$method = gsub('.*-(.*)-.*-.*-.*-prcurve.tsv', '\\1', basename(ff))
    df$sample = gsub('.*-.*-(.*)-.*-.*-prcurve.tsv', '\\1', basename(ff))
    df$region = gsub('.*-.*-.*-(.*)-.*-prcurve.tsv', '\\1', basename(ff))
    df$eval = gsub('.*-.*-.*-.*-(.*)-prcurve.tsv', '\\1', basename(ff))
    df$min.cov = .5
    df
  })
  pr.df = do.call(rbind, df)
  utils::write.table(pr.df, file=paste0(cmd.args$out$value, '-prcurve.tsv'), quote=FALSE, sep='\t', row.names=FALSE)

  ## Combined figure of the best F1
  ## library(sveval)
  ## library(ggplot2)
  ## library(dplyr)
  ## library(RColorBrewer)

  ## Relabels columns and set orders
  ##   reg.labels contains labels to associate with region names,
  ##     also defining the order to show in the graphs
  relabel <- function(df, reg.labels=c(all='all',
                                       rep='repeat',
                                       nonrep='non-repeat',
                                       hc='high-confidence',
                                       called='called in SMRT-SV v2',
                                       nocalls='not called in SMRT-SV v2',
                                       conf='high-confidence')){
    if('method' %in% colnames(df)){
      df$method = factor(df$method, levels=unique(df$method))
    }
    ## Types
    if('type' %in% colnames(df)){
      df$type = factor(df$type, levels=c('Total', 'INS', 'DUP', 'DEL', 'INV', 'BND'))
    }
    ## Region
    if('region' %in% colnames(df)){
      regions = unique(df$region)
      for(missing.label in setdiff(regions, names(reg.labels))){
        reg.labels[[missing.label]] = missing.label
      }
      df$region = factor(df$region, levels=names(reg.labels), labels=reg.labels)
    }
    ## Evaluation metric
    if('eval' %in% colnames(df)){
      df$eval=factor(df$eval, levels=c('call','geno'),
                     labels=c('presence', 'genotype'))
    }
    ## Sizes
    if('size' %in% colnames(df)){
      sizes = unique(df$size)
      sizes = sizes[order(as.numeric(gsub('.*,(.*)]', '\\1', sizes)))]
      sizes.l = gsub('\\((.*),Inf]', '>\\1', sizes)
      sizes.l = gsub('e\\+03', 'K', sizes.l)
      sizes.l = gsub('e\\+04', '0K', sizes.l)
      sizes.l = gsub('e\\+05', '00K', sizes.l)
      df$size = factor(df$size, levels=sizes, labels=sizes.l)
    }
    return(df)
  }

  ## Remove Total, order and relabel
  pr.df = pr.df[which(pr.df$type != 'Total'),]
  pr.df = pr.df[order(as.numeric(pr.df$qual)),]
  pr.df = relabel(pr.df)


  ## If there are no SV of a certain type, remove it from the graph
  for(svtype in c('DEL', 'INS', 'DUP', 'INV', 'BND')){
    if(!any(pr.df$type==svtype & !is.na(pr.df$F1))) {
      pr.df = pr.df[which(pr.df$type!=svtype),]
    }
  }
  
  ## Sum the TP/FP/FN across samples and recompute precision/recall/F1
  pr.df = pr.df %>% dplyr::group_by(.data$exp, .data$type, .data$qual,
                                    .data$method, .data$region, .data$eval) %>%
    dplyr::select(.data$TP, .data$TP.baseline, .data$FN, .data$FP) %>%
    dplyr::summarize_all(sum) %>%
    prf

  ## est F1 for each method. for the bar plots
  eval.f1 = pr.df %>%
    dplyr::group_by(.data$exp, .data$method, .data$type, .data$region, .data$eval) %>%
    dplyr::arrange(dplyr::desc(.data$F1)) %>%
    dplyr::do(utils::head(.data, 1))

  ggp = list()

  ## F1 score bar graph
  recall = precision = F1 = type = region = eval = method = NULL
  ggp$f1 = eval.f1 %>% 
    ggplot(aes(x=region, y=F1, fill=method, alpha=eval, group=method)) +
    geom_bar(stat='identity', color='black', position=position_dodge()) +
    facet_grid(exp~type, scales='free', space='free') +
    scale_fill_brewer(palette='Set1') + 
    scale_alpha_manual(name='SV evaluation', values=c(.5,1)) + 
    theme_bw() + ylim(0,1) + 
    labs(x='Genomic regions', y='Best F1', fill='Method')

  ## precision score bar graph
  ggp$prec = eval.f1 %>% 
    ggplot(aes(x=region, y=precision, fill=method, alpha=eval, group=method)) +
    geom_bar(stat='identity', color='black', position=position_dodge()) +
    facet_grid(exp~type, scales='free', space='free') +
    scale_fill_brewer(palette='Set1') + 
    scale_alpha_manual(name='SV evaluation', values=c(.5,1)) + 
    theme_bw() + ylim(0,1) + 
    labs(x='Genomic regions', y='precision', fill='Method')

  ## recall score bar graph
  ggp$recall = eval.f1 %>% 
    ggplot(aes(x=region, y=recall, fill=method, alpha=eval, group=method)) +
    geom_bar(stat='identity', color='black', position=position_dodge()) +
    facet_grid(exp~type, scales='free', space='free') +
    scale_fill_brewer(palette='Set1') + 
    scale_alpha_manual(name='SV evaluation', values=c(.5,1)) + 
    theme_bw() + ylim(0,1) + 
    labs(x='Genomic regions', y='recall', fill='Method')

  ## precision-recall curves for each experiment 'exp'
  for(exp.c in unique(pr.df$exp)){
    ggp$pr = pr.df %>% dplyr::filter(.data$exp==exp.c) %>% dplyr::arrange(.data$qual) %>% 
      ggplot(aes(x=recall, y=precision, color=method, linetype=region)) +
      geom_point(aes(shape=region)) + geom_path() + theme_bw() +
      xlim(0,1) + ylim(0,1) +
      ggtitle(exp.c) + 
      facet_grid(eval~type)
    names(ggp)[length(ggp)] = paste0('pr_', exp.c)
  }

  ## save graphs to PDF
  grDevices::pdf(paste0(cmd.args$out$value, '.pdf'), 8, 6)
  tmp = lapply(ggp, print)
  grDevices::dev.off()

}
