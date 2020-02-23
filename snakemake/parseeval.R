## Extract some info from the RData file and write TSV files
## Arguments
## 1: R file with the output of svevalOl as a 'eval.o' object (computed in sveval.R wrapper)
## 2: output PR curve TSV file
## 3: output per-size TSV file
args = commandArgs(TRUE)

library(sveval)

## load RData object
load(args[1])

## PR curve
write.table(eval.o$curve, file=args[2], sep='\t', quote=FALSE, row.names=FALSE)

## per-size estimates
df = plot_persize(eval.o, plot=FALSE,
                  size.breaks=c(50,100,200,300,400,600,800,
                                1000,2500,5000,Inf))
write.table(df, file=args[3], sep='\t', quote=FALSE, row.names=FALSE)
