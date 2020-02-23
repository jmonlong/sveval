## wrapper to read TSV files, merge them and make some graphs.
args = commandArgs(TRUE)
## arg 1: file name for the output PDF file 
out.pdf = args[1]
## arg 2: the list of TSV files to read and merge
all.tsvs = unlist(strsplit(args[2], ' '))

## prefix of the output file. used for the merged TSV file
out.prefix = gsub('.pdf$', '', basename(out.pdf))

## Merge persize files
tsvs = grep('persize', all.tsvs, value=TRUE)
df = lapply(tsvs, function(ff){
  df = read.table(ff, as.is=TRUE, header=TRUE)
  ## extract information from the filename
  df$exp = gsub('tsv/(.*)-.*-.*-.*-.*-persize.tsv', '\\1', ff)
  df$method = gsub('tsv/.*-(.*)-.*-.*-.*-persize.tsv', '\\1', ff)
  df$sample = gsub('tsv/.*-.*-(.*)-.*-.*-persize.tsv', '\\1', ff)
  df$region = gsub('tsv/.*-.*-.*-(.*)-.*-persize.tsv', '\\1', ff)
  df$eval = gsub('tsv/.*-.*-.*-.*-(.*)-persize.tsv', '\\1', ff)
  df$min.cov = .5
  df
})
df = do.call(rbind, df)
write.table(df, file=paste0('tsv/', out.prefix, '-persize.tsv'), quote=FALSE, sep='\t', row.names=FALSE)

## Merge PR curve files
tsvs = grep('prcurve', all.tsvs, value=TRUE)
df = lapply(tsvs, function(ff){
  df = read.table(ff, as.is=TRUE, header=TRUE)
  ## extract information from the filename
  df$exp = gsub('tsv/(.*)-.*-.*-.*-.*-prcurve.tsv', '\\1', ff)
  df$method = gsub('tsv/.*-(.*)-.*-.*-.*-prcurve.tsv', '\\1', ff)
  df$sample = gsub('tsv/.*-.*-(.*)-.*-.*-prcurve.tsv', '\\1', ff)
  df$region = gsub('tsv/.*-.*-.*-(.*)-.*-prcurve.tsv', '\\1', ff)
  df$eval = gsub('tsv/.*-.*-.*-.*-(.*)-prcurve.tsv', '\\1', ff)
  df$min.cov = .5
  df
})
pr.df = do.call(rbind, df)
write.table(pr.df, file=paste0('tsv/', out.prefix, '-prcurve.tsv'), quote=FALSE, sep='\t', row.names=FALSE)


## Combined figure of the best F1
library(sveval)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

## Relabels columns and set orders
relabel <- function(df, nonrep=c('nonrep', 'hc')){
  if('method' %in% colnames(df)){
    df$method = factor(df$method, levels=unique(df$method))
  }
  ## Types
  if('type' %in% colnames(df)){
    df$type = factor(df$type, levels=c('Total', 'INS', 'DEL', 'INV'))
  }
  ## Region
  if('region' %in% colnames(df)){
    reg.l = c('all','repeat', 'non-repeat', 'called in SMRT-SV v2',
              'not called in SMRT-SV v2')
    if(nonrep[1] == 'nonrep'){
      reg.l[3] = 'non-repeat'
    } else if(nonrep[1] == 'hc'){
      reg.l[3] = 'high-confidence'
    }
    df$region = factor(df$region, levels=c('all','rep', 'nonrep', 'called', 'nocalls'),
                       labels=reg.l)
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

## Read evaluation results
pr.df = pr.df %>% filter(type!='Total') %>% arrange(qual)
pr.df = relabel(pr.df)

## If there are no inversion, remove this type from the graph
if(!any(pr.df$type=='INV' & !is.na(pr.df$F1))) {
  pr.df = subset(pr.df, type!='INV')
}

## Sum the TP/FP/FN across samples and recompute precision/recall/F1
pr.df = pr.df %>% group_by(exp, exp, type, qual, method, region, eval) %>%
  select(TP, TP.baseline, FN, FP) %>% summarize_all(sum)
pr.df = prf(pr.df)

## est F1 for each method. for the bar plots
eval.f1 = pr.df %>% group_by(exp, exp, method, type, region, eval) %>%
  arrange(desc(F1)) %>% do(head(., 1))

ggp = list()

## F1 score bar graph
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
  ggp$pr = pr.df %>% filter(exp==exp.c) %>% arrange(qual) %>% 
    ggplot(aes(x=recall, y=precision, color=method, linetype=region)) +
    geom_point(aes(shape=region)) + geom_path() + theme_bw() +
    xlim(0,1) + ylim(0,1) +
    ggtitle(exp.c) + 
    facet_grid(eval~type)
  names(ggp)[length(ggp)] = paste0('pr_', exp.c)
}

## save graphs to PDF
pdf(out.pdf, 8, 6)
tmp = lapply(ggp, print)
dev.off()

