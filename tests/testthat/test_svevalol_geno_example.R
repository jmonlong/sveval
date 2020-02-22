context('Genotype evaluation in some specific situations')

test_that("Example VCF", {
  res = svevalOl('../example-calls.vcf', '../example-truth.vcf', geno.eval=TRUE,
                 min.size=2, merge.hets=TRUE, stitch.hets=TRUE, method='bipartite')
  expect_true(all(res$TP>0))
  expect_true(all(res$TP.baseline>0))
  expect_true(all(res$FP==0))
  expect_true(all(res$FN==0))
})
