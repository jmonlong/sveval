context('Overlap function')
library(GenomicRanges)

test_that("REF/ALT VCFs", {
  calls = readSVvcf('../calls.s0.vcf', check.inv=TRUE)
  truth = readSVvcf('../truth.refalt.vcf', check.inv=TRUE)
  res = svOverlap(truth, calls)
  expect_gt(length(res$queryHits), 0)
  expect_gt(sum(res$queryOl), 0)
  expect_gt(length(res$subjectHits), 0)
  expect_gt(sum(res$subjectOl), 0)
})

test_that("Overlap simple insertions extended with simple repeat annotation", {
  sv1.gr = GRanges('x', IRanges(10, width=1), size=100, type='INS')
  sv2.gr = GRanges('x', IRanges(100, width=1), size=100, type='INS')
  res = svOverlap(sv1.gr, sv2.gr)
  expect_true(is.null(res))
  simprep.gr = GRanges('x', IRanges(10, width=100))
  res = svOverlap(sv1.gr, sv2.gr, simprep=simprep.gr)
  expect_true(!is.null(res))
})


