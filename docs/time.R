library(devtools)
load_all()

hgsvc.truth = '~/Documents/ucsc/sv/HGSVC.haps.vcf.gz'
hgsvc.calls = "~/Documents/ucsc/sv/mergedcatalog/sveval/hgsvc-hsvlr_srdedup17_aug_18a5600_gssw_gaf-HG00514.vcf.gz"

'VCF reading'
system.time((truth.gr = readSVvcf(hgsvc.truth)))
system.time((calls.gr = readSVvcf(hgsvc.calls, qual.field='QUAL')))

'Calling evaluation'
system.time((eval.o = svevalOl(calls.gr, truth.gr)))

'Genotype evaluation (coverage)'
system.time((eval.o = svevalOl(calls.gr, truth.gr, geno.eval=TRUE, stitch.hets=TRUE, merge.hets=TRUE)))

'Genotype evaluation (bipartite)'
system.time((eval.o = svevalOl(calls.gr, truth.gr, geno.eval=TRUE, stitch.hets=TRUE, merge.hets=TRUE, method='bipartite')))

'VCF reading of SNV/indels/SV'
pepper = "~/Documents/ucsc/sv/sveval/HG002_R941_GIAB_minION_guppy324_PEPPER_HP_candidates_merged.fixed.vcf.gz"
system.time((pepper.gr = readSVvcf(pepper)))


## ## Debug differences
## library(GenomicRanges)

## truth.rcpp = readSVvcf(hgsvc.truth)
## calls.rcpp = readSVvcf(hgsvc.calls, qual.field='QUAL')
## eval.rcpp = svevalOl(calls.rcpp.f, truth.rcpp)

## truth.master = readSVvcf(hgsvc.truth)
## calls.master = readSVvcf(hgsvc.calls, right.trim=FALSE, qual.field='QUAL')
## eval.master = svevalOl(calls.master.f, truth.master)

## eval.rcpp$eval
## eval.master$eval

## truth.rcpp.f = subset(truth.rcpp, ac>0)
## length(truth.rcpp.f)
## length(truth.master)

## calls.rcpp.f = subset(calls.rcpp, ac>0)
## length(calls.rcpp.f)
## calls.master.f = subset(calls.master, size>=10)
## length(calls.master.f)

## table(calls.rcpp.f$type)
## table(calls.master.f$type)

## calls.rcpp.f$svid = paste(calls.rcpp.f$type, as.character(seqnames(calls.rcpp.f)), start(calls.rcpp.f))
## calls.master.f$svid = paste(calls.master.f$type, as.character(seqnames(calls.master.f)), start(calls.master.f))

## calls.rcpp.f$svid = paste(calls.rcpp.f$type, as.character(calls.rcpp.f))
## calls.master.f$svid = paste(calls.master.f$type, as.character(calls.master.f))


## svids.rcpp = setdiff(calls.rcpp.f$svid, calls.master.f$svid)
## svids.master = setdiff(calls.master.f$svid, calls.rcpp.f$svid)

## all(calls.rcpp.f$svid == calls.master.f$svid)
## all(calls.rcpp.f$size == calls.master.f$size)
## all(calls.rcpp.f$type == calls.master.f$type)
## all(calls.rcpp.f$qual == calls.master.f$QUAL)
## all(ifelse(calls.rcpp.f$ac==1, 'het', 'hom') == calls.master.f$GT)

## ex = sample(svids.rcpp, 1)
## subset(calls.rcpp.f, svid == ex)
## subsetByOverlaps(calls.master.f, subset(calls.rcpp.f, svid == ex))

## subset(calls.master.f, type=="INS" & width(calls.master.f)>3)

## calls.rcpp = readSVvcf(hgsvc.calls, qual.field='QUAL', out.fmt='df')


## (inserted) sequence alignment
library(GenomicRanges)
svs = read.table('~/Documents/ucsc/topmed/mesa_batch12/svs.batch1_2.seq.tsv.gz', header=TRUE, as.is=TRUE)
svs = makeGRangesFromDataFrame(svs, keep.extra.columns = TRUE)
## quick clustering to split dataset for analysis
cl.gr = reduce(svs, min.gapwidth=100)
cl.gr$n = countOverlaps(cl.gr, svs)
cl.gr$n.large.ins = countOverlaps(cl.gr, subset(svs, type=='INS' & size > 2000))
ex.gr = subset(cl.gr, n>10 & n.large.ins>10)[1]
svs.ex = subsetByOverlaps(svs, ex.gr)
svs.ex = sample(svs.ex, 30)

summary(svs.ex$size)

system.time((ol = svOverlap(svs.ex, svs.ex, ins.seq.comp=TRUE)))

## R-bioc
##   user  system elapsed 
## 78.060   0.016  78.140 

## edlib
##  user  system elapsed 
## 0.372   0.000   0.372 

## 200x faster on ~3Kbp insertions
