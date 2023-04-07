context('SV clustering')

test_that("Calls and truth calls cluster together", {
  calls.gr = readSVvcf('../calls.s0.vcf')
  truth.gr = readSVvcf('../truth.refalt.vcf')
  cl.gr = clusterSVs(c(calls.gr, truth.gr))
  expect_true(any(table(cl.gr$svsite)>1))
})
