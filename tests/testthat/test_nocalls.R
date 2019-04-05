context('No-calls identification')

test_that("Input with symbolic VCF representation", {
  nocalls.idx = findNocalls('../calls.s0.vcf', '../truth.symb.vcf')
  expect_gt(nrow(nocalls.idx), 0)
  expect_gt(ncol(nocalls.idx), 0)
})

test_that("No no-calls variants", {
  nocalls.idx = findNocalls('../empty.vcf', '../truth.symb.vcf')
  expect_equal(nrow(nocalls.idx), 0)
})
