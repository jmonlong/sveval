context('Overlap function')

test_that("REF/ALT VCFs", {
  calls = readSVvcf('../calls.s0.vcf', check.inv=TRUE)
  truth = readSVvcf('../truth.refalt.vcf', check.inv=TRUE)
  res = svOverlap(truth, calls)
  expect_gt(length(res$queryHits), 0)
  expect_gt(sum(res$queryOl), 0)
  expect_gt(length(res$subjectHits), 0)
  expect_gt(sum(res$subjectOl), 0)
})

