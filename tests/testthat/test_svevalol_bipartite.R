context('Overlap-based comparison for genotype evaluation using bipartite clustering to match variants')

test_that("ALT/REF inputs and output in file", {
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', outfile='temp.tsv', geno.eval=TRUE, method='bipartite')
  res = read.table('temp.tsv', header=TRUE, as.is=TRUE, sep='\t')
  expect_gt(nrow(res), 0)
  expect_true(any(res$TP>0))
  expect_true(any(res$TP.baseline>0))
  file.remove('temp.tsv')
})

test_that("Stitch and merge hets", {
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', geno.eval=TRUE, merge.hets=TRUE, stitch.hets=TRUE, ins.seq.comp=TRUE, method='bipartite')
  res = res$eval
  expect_gt(nrow(res), 0)
  expect_true(any(res$TP>0))
  expect_true(any(res$TP.baseline>0))
})

test_that("Output BED files etc", {
  res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf', out.bed.prefix='tempfortest', geno.eval=TRUE, method='bipartite')
  res = res$eval
  expect_gt(nrow(res), 0)
  expect_true(any(res$TP>0))
  expect_true(any(res$TP.baseline>0))
  file.remove(list.files('.', 'tempfortest'))
})

test_that("Input with symbolic VCF representation", {
  res = svevalOl('../calls.s0.vcf', '../truth.symb.vcf', geno.eval=TRUE, method='bipartite')
  expect_gt(nrow(res$eval), 0)
  expect_true(any(as.matrix(res$eval[,2:5])>0))
})

test_that("Empty inputs", {
  calls.gr = readSVvcf('../calls.s0.vcf')
  truth.gr = readSVvcf('../truth.refalt.vcf')
  ## Empty calls
  res = svevalOl('../empty.vcf', truth.gr, min.size=20, geno.eval=TRUE, method='bipartite')
  expect_true(all(res$eval$recall == 0))
  expect_true(all(res$curve$recall == 0))
  ## Empty truth set
  expect_error(svevalOl(calls.gr, '../empty.vcf', geno.eval=TRUE, method='bipartite'), "no SVs")
  ## One type missing
  res = svevalOl(calls.gr[which(calls.gr$type=='DEL')], truth.gr, min.size=20, geno.eval=TRUE, method='bipartite')
  res = res$eval
  expect_true(any(as.matrix(res[,2:5])>0))
})
