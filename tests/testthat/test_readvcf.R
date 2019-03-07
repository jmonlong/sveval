context('Read VCFs containing SVs')

test_that("ALT/REF VCF", {
  vcf = readSVvcf('../truth.refalt.vcf')
  expect_gt(length(vcf), 0)
})

test_that("VCF with symbolic SVs", {
  vcf = readSVvcf('../truth.symb.vcf')
  expect_gt(length(vcf), 0)
})

test_that("VCF with no QUAL but GQ field", {
  vcf = readSVvcf('../delly.vcf')
  expect_gt(length(vcf), 0)
  expect_gt(sum(vcf$QUAL), 0)
})

test_that("VCF with existing QUAL and another FORMAT field", {
  vcf = readSVvcf('../calls.s0.vcf')
  vcf2 = readSVvcf('../calls.s0.vcf', qual.field='DP')
  expect_true(any(vcf$QUAL != vcf2$QUAL))
})

