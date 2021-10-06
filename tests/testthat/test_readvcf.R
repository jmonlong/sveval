context('Read VCFs containing SVs')

test_that("ALT/REF VCF", {
  vcf = readSVvcf('../truth.refalt.vcf')
  expect_gt(length(vcf), 0)
})

test_that("ALT/REF VCF checking for inversions", {
  vcf = readSVvcf('../truth.refalt.vcf', check.inv=TRUE)
  expect_gt(length(vcf), 0)
})

test_that("VCF with symbolic SVs and different samples", {
  vcf = readSVvcf('../truth.symb.vcf')
  expect_gt(length(vcf), 0)
  vcf2 = readSVvcf('../truth.symb.vcf', sample.name='s1')
  expect_gt(length(vcf2), 0)
  vcf = subset(vcf, ac>0)
  vcf2 = subset(vcf2, ac>0)
  expect_true(length(vcf) != length(vcf2))
})

test_that("VCF gzipped", {
  vcf = readSVvcf('../delly.vcf.gz')
  expect_gt(length(vcf), 0)
  expect_gt(sum(vcf$qual), 0)
  
})

test_that("VCF bgzipped", {
  vcf = readSVvcf('../delly.vcf.bgz')
  expect_gt(length(vcf), 0)
  expect_gt(sum(vcf$qual), 0)
  expect_gt(sum(vcf$filter == 'PASS'), 0)
})

test_that("VCF with no QUAL but GQ field", {
  vcf = readSVvcf('../delly.vcf')
  expect_gt(length(vcf), 0)
  expect_gt(sum(vcf$qual), 0)
})

test_that("VCF with existing QUAL and another FORMAT field", {
  vcf = readSVvcf('../calls.s0.vcf')
  vcf2 = readSVvcf('../calls.s0.vcf', qual.field='DP')
  expect_true(any(vcf$qual != vcf2$qual))
})

test_that("VCF with missing quality filter", {
  vcf = readSVvcf('../calls.s0.vcf')
  vcf2 = readSVvcf('../calls.s0.vcf', qual.field='QAUL')
  expect_true(any(vcf$qual == vcf2$qual))
})

test_that("ALT/REF VCF with VCF output", {
  vcf = readSVvcf('../truth.refalt.vcf', out.fmt='vcf')
  expect_gt(length(vcf), 0)
})

test_that("VCF with symbolic SVs with VCF output", {
  vcf = readSVvcf('../truth.symb.vcf', out.fmt='vcf')
  expect_gt(length(vcf), 0)
})

test_that("ALT/REF VCF with data.frame output", {
  vcf = readSVvcf('../truth.refalt.vcf', out.fmt='df')
  expect_gt(length(vcf), 0)
})

test_that("BND and TRA types are retrieved", {
  vcf = readSVvcf('../truth.symb.vcf', out.fmt='df', other.field='CHR2')
  expect_gt(length(vcf), 0)
  expect_true(any(vcf$type=='BND'))
  expect_true(any(vcf$type=='TRA'))
  expect_true(any(!is.na(vcf$CHR2)))
  expect_true(any(!is.na(vcf$end2)))
})

test_that("Saving specified INFO field", {
  vcf = readSVvcf('../truth.symb.vcf', out.fmt='df', other.field=c('AC','FREQ'))
  expect_gt(length(vcf), 0)
})

