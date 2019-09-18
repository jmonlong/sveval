context('Frequency annotation')

test_that("Reading VCFs", {
  res = freqAnnotate('../calls.s0.vcf', '../truth.symb.vcf', freq.field='FREQ', out.freq.field='FREQ')
  expect_true(mean(VariantAnnotation::info(res)$FREQ)>0)
})

test_that("Readinf VCFs and writing output VCF", {
  res = freqAnnotate('../calls.s0.vcf', '../truth.symb.vcf', freq.field='FREQ',
                     out.vcf='temp.vcf')
  expect_true(file.remove('temp.vcf'))
})


