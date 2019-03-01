context('Overlap-based comparison')

test_that("ALT/REF inputs and output in file", {
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', outfile='temp.tsv')
  res = read.table('temp.tsv', header=TRUE, as.is=TRUE, sep='\t')
  expect_gt(nrow(res), 0)
  expect_true(all(res$TP>0))
  expect_true(all(res$TP.baseline>0))
  file.remove('temp.tsv')
})

test_that("Output BED files etc", {
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', out.bed.prefix='tempfortest')
  res = res$eval
  expect_gt(nrow(res), 0)
  expect_true(all(res$TP>0))
  expect_true(all(res$TP.baseline>0))
  file.remove(list.files('.', 'tempfortest'))
})

test_that("Input with symbolic VCF representation", {
  res = svevalOl('../calls.s0.vcf', '../truth.symb.vcf')
  expect_gt(nrow(res$eval), 0)
  expect_true(any(as.matrix(res$eval[,2:5])>0))
  pdf('temp.pdf')
  print(plot_prcurve(res$curve))
  dev.off()
  expect_true(file.remove('temp.pdf'))
})

test_that("Filters", {
  res.all = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf',  min.size=0)
  res.all = res.all$eval
  ## BED file
  bed = data.frame(chr='x', start=c(1e5, 7e5), end=c(5e5, 1e6))
  write.table(bed, file='temp.bed', row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', bed.regions='temp.bed', min.size=0)
  res = res$eval
  expect_gt(nrow(res), 0)
  expect_true(all(res$TP>0))
  expect_true(all(res$TP.baseline>0))
  expect_true(all(as.matrix(res[1:3,2:5])<=as.matrix(res.all[1:3,2:5])))
  file.remove('temp.bed')
  ## Small variants
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', min.size=20)
  res = res$eval
  expect_gt(nrow(res), 0)
  expect_true(all(res$TP>0))
  expect_true(all(res$TP.baseline>0))
  expect_true(all(as.matrix(res[,2:5])<=as.matrix(res.all[,2:5])))
})


test_that("Sequence comparison for insertions", {
  calls.gr = readSVvcf('../calls.s0.vcf', keep.ins.seq=TRUE)
  ## Subset insertions to speed up the test
  ins.idx = sample(which(calls.gr$type=='INS'), 10)
  del.idx = which(calls.gr$type=='DEL')
  calls.gr = calls.gr[c(ins.idx, del.idx)]
  ## Run evaluation
  res = svevalOl(calls.gr, '../truth.refalt.vcf', min.size=20, ins.seq.comp=TRUE)
  res = res$eval
  expect_gt(nrow(res), 0)
  expect_true(any(as.matrix(res[1:3,2:5])>0))
  expect_true(all(res$TP.baseline[1:3]>0))
})


test_that("Empty inputs", {
  calls.gr = readSVvcf('../calls.s0.vcf')
  truth.gr = readSVvcf('../truth.refalt.vcf')
  ## Empty calls
  res = svevalOl('../empty.vcf', truth.gr, min.size=20)
  res = res$eval
  expect_true(all(is.na(as.matrix(res[,2:5]))))
  ## Empty truth set
  expect_error(svevalOl(calls.gr, '../empty.vcf'), "no SVs")
  ## One type missing
  res = svevalOl(calls.gr[which(calls.gr$type=='DEL')], truth.gr, min.size=20)
  res = res$eval
  expect_true(any(as.matrix(res[,2:5])>0))
})
