context('Frequency annotation')

test_that("Reading VCFs with inv and dup", {
  res = freqAnnotate('../calls.s0.vcf', '../truth.symb.vcf', freq.field='FREQ', out.freq.field='FREQ', check.inv=TRUE)
  expect_true(mean(VariantAnnotation::info(res)$FREQ)>0)
})

test_that("Error when missing frequency field", {
  expect_error({freqAnnotate('../calls.s0.vcf', '../truth.refalt.vcf', freq.field='FREQ')}, 'not an INFO field')
})

test_that("Readinf VCFs and writing output VCF", {
  res = freqAnnotate('../calls.s0.vcf', '../truth.symb.vcf', freq.field='FREQ',
                     out.vcf='temp.vcf')
  expect_true(file.remove('temp.vcf'))
})


