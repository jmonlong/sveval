context('Overlap function')

test_that("REF/ALT VCFs", {
  calls = readSVvcf('../calls.s0.vcf', check.inv=TRUE)
  truth = readSVvcf('../truth.refalt.vcf', check.inv=TRUE)
  res = svOverlap(truth, calls)
  expect_gt(length(res$query), 0)
  expect_gt(sum(res$query$cov.prop), 0)
  expect_gt(length(res$subject), 0)
  expect_gt(sum(res$subject$cov.prop), 0)
})

