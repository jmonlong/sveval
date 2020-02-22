context('Genotype evaluation in some specific situations')

test_that("Example VCF", {
  res = svevalOl('../example-calls.vcf', '../example-truth.vcf', geno.eval=TRUE,
                 min.size=2, merge.hets=FALSE, stitch.hets=FALSE, method='bipartite')
  ## multi-allele correctly matched
  expect_true('x:500-507' %in% as.character(res$svs$DEL$TP.baseline))
  ## similar deletion correctly priorized
  expect_true('x:1000-1062' %in% as.character(res$svs$DEL$TP.baseline))
})
