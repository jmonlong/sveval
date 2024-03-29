context('Overlap-based comparison')
library(GenomicRanges)

test_that("ALT/REF inputs and output in file", {
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', outfile='temp.tsv', out.bed.prefix='tempfortest')
  res = read.table('temp.tsv', header=TRUE, as.is=TRUE, sep='\t')
  expect_gt(nrow(res), 0)
  expect_gt(sum(res$TP>0), 2)
  expect_gt(sum(res$TP.baseline>0), 2)
  file.remove('temp.tsv')
})

test_that("Input with symbolic VCF representation", {
  simprep.gr = GRanges('x', IRanges(10, width=100))
  res = svevalOl('../calls.s0.vcf', '../truth.symb.vcf', simprep=simprep.gr)
  expect_gt(nrow(res$eval), 0)
  expect_gt(sum(res$eval$TP>0), 2)
  expect_gt(sum(res$eval$TP.baseline>0), 2)
  expect_gt(sum(res$eval$TP.baseline>0), 2)
})

test_that("Filters", {
  res.all = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf',  min.size=0)
  res.all = res.all$eval
  ## BED file
  bed = data.frame(chr='x', start=c(3000, 1480000), end=c(9000, 1490000))
  write.table(bed, file='temp.bed', row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', bed.regions='temp.bed', min.size=0)
  res = res$eval
  expect_gt(nrow(res), 0)
  expect_gt(sum(res$TP>0), 2)
  expect_gt(sum(res$TP.baseline>0), 2)
  res.m = merge(res[,c('type', 'TP')], res.all[,c('type', 'TP')], by='type')
  expect_true(all(res.m$TP.x<=res.m$TP.y))
  file.remove('temp.bed')
  ## Small variants
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', min.size=20)
  res = res$eval
  expect_gt(nrow(res), 0)
  expect_gt(sum(res$TP>0), 2)
  expect_gt(sum(res$TP.baseline>0), 2)
  res.m = merge(res[,c('type', 'TP')], res.all[,c('type', 'TP')], by='type')
  expect_true(all(res.m$TP.x<=res.m$TP.y))
})

test_that("Sequence comparison for insertions and deletions", {
  calls.gr = readSVvcf('../calls.s0.vcf', keep.ins.seq=TRUE, keep.ref.seq=TRUE)
  ## Subset variants to speed up the test
  ins.idx = sample(which(calls.gr$type=='INS' & calls.gr$ac>0), 10)
  del.idx = sample(which(calls.gr$type=='DEL' & calls.gr$ac>0), 10)
  calls.gr = calls.gr[c(ins.idx, del.idx)]
  ## Run evaluation
  res = svevalOl(calls.gr, '../truth.refalt.vcf', min.size=20, ins.seq.comp=TRUE, range.seq.comp=TRUE)
  res = res$eval
  expect_gt(nrow(res), 0)
  expect_true(res$TP[which(res$type=='INS')]>0)
  expect_true(res$TP[which(res$type=='DEL')]>0)
})

test_that("Sequence comparison error when no seq available", {
  calls.gr = readSVvcf('../calls.s0.vcf', )
  ## Subset variants to speed up the test
  ins.idx = sample(which(calls.gr$type=='INS' & calls.gr$ac>0), 10)
  del.idx = sample(which(calls.gr$type=='DEL' & calls.gr$ac>0), 10)
  calls.gr = calls.gr[c(ins.idx, del.idx)]
  ## Run evaluation
  expect_error(svevalOl(calls.gr, '../truth.refalt.vcf', min.size=20, ins.seq.comp=TRUE, range.seq.comp=TRUE), 'Missing sequence')
})

test_that("Sequence comparison for insertions and deletions (multicore)", {
  calls.gr = readSVvcf('../calls.s0.vcf', keep.ins.seq=TRUE, keep.ref.seq=TRUE)
  ## Subset insertions to speed up the test
  ins.idx = sample(which(calls.gr$type=='INS' & calls.gr$ac>0), 10)
  del.idx = sample(which(calls.gr$type=='DEL' & calls.gr$ac>0), 10)
  calls.gr = calls.gr[c(ins.idx, del.idx)]
  ## Run evaluation
  res = svevalOl(calls.gr, '../truth.refalt.vcf', min.size=20, ins.seq.comp=TRUE, range.seq.comp=TRUE, nb.cores=2)
  res = res$eval
  expect_gt(nrow(res), 0)
  expect_true(res$TP[which(res$type=='INS')]>0)
  expect_true(res$TP[which(res$type=='DEL')]>0)
})


test_that("Empty inputs", {
  calls.gr = readSVvcf('../calls.s0.vcf')
  truth.gr = readSVvcf('../truth.refalt.vcf')
  ## Subset variants to speed up the test
  ins.idx = sample(which(calls.gr$type=='INS' & calls.gr$ac>0), 10)
  del.idx = sample(which(calls.gr$type=='DEL' & calls.gr$ac>0), 10)
  calls.gr = calls.gr[c(ins.idx, del.idx)]
  ## Empty calls
  expect_warning((res = svevalOl('../empty.vcf', truth.gr, min.size=20)), "no SVs")
  expect_true(all(res$eval$recall == 0))
  expect_true(all(res$curve$recall == 0))
  ## Empty truth set
  expect_error(svevalOl(calls.gr, '../empty.vcf'), "no SVs")
  ## One type missing
  res = svevalOl(calls.gr[which(calls.gr$type=='DEL')], truth.gr, min.size=20)
  res = res$eval
  expect_true(any(res$TP==0))
  expect_true(any(res$TP>0))
  ## Calls are all homozygous ref
  calls.gr$ac = 0
  expect_warning(svevalOl(calls.gr, truth.gr), "no SVs")
})
