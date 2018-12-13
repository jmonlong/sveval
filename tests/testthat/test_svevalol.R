context('Overlap-based comparison')

test_that("ALT/REF inputs and output in file", {
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', outfile='temp.tsv')
  res = read.table('temp.tsv', header=TRUE, as.is=TRUE, sep='\t')
  expect_gt(nrow(res), 0)
  expect_true(all(res$TP>0))
  expect_true(all(res$TP.baseline>0))
  file.remove('temp.tsv')
})

test_that("Input with symbolic VCF representation", {
  res = svevalOl('../calls.s0.vcf', '../truth.symb.vcf')
  expect_gt(nrow(res), 0)
  expect_true(all(res$TP>0))
  expect_true(all(res$TP.baseline>0))
})

test_that("Filters", {
  res.all = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf')
  ## BED file
  bed = data.frame(chr='x', start=c(1e5, 7e5), end=c(5e5, 1e6))
  write.table(bed, file='temp.bed', row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', bed.regions='temp.bed')
  expect_gt(nrow(res), 0)
  expect_true(all(res$TP>0))
  expect_true(all(res$TP.baseline>0))
  expect_true(all(as.matrix(res[,2:5])<=as.matrix(res.all[,2:5])))
  file.remove('temp.bed')
  ## Small variants
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', min.size=20)
  expect_gt(nrow(res), 0)
  expect_true(all(res$TP>0))
  expect_true(all(res$TP.baseline>0))
  expect_true(all(as.matrix(res[,2:5])<=as.matrix(res.all[,2:5])))
})


test_that("Sequence compatison for insertions", {
  calls.gr = readSVvcf('../calls.s0.vcf', keep.ins.seq=TRUE)
  ## Subset insertions to speed up the test
  ins.idx = sample(which(calls.gr$type=='INS'), 10)
  del.idx = which(calls.gr$type=='INS')
  calls.gr = calls.gr[c(ins.idx, del.idx)]
  ## Run evaluation
  res = svevalOl(calls.gr, '../truth.refalt.vcf', min.size=20, ins.seq.comp=TRUE)
  expect_gt(nrow(res), 0)
  expect_true(all(res$TP>0))
  expect_true(all(res$TP.baseline>0))
})

