context('Read VCFs containing SVs')

test_that("Multisamples with ALT/REF", {
  vcf = readSVvcf.multisamps('../truth.refalt.vcf')
  expect_gt(length(vcf), 0)
  expect_true(any(vcf$af>0))
  vcf = readSVvcf.multisamps('../truth.refalt.vcf', keep.ins.seq=TRUE, keep.ref.seq=TRUE,
                             check.inv=TRUE, keep.ids=TRUE, out.fmt='df')
  expect_gt(length(vcf), 0)
  expect_true(any(vcf$af>0))
})

test_that("Multisamples with ALT/REF", {
  df = readSVvcf.multisamps('../truth.refalt.vcf', out.fmt='df', keep.ids=TRUE, min.sv.size=100)
  sites = list(df$svid[1], df$svid[2], df$svid[3:10], df$svid[11:40])
  ac = countAlleles('../truth.refalt.vcf', sv.sites=sites)
  ac.tot = sapply(sites, function(svids) sum(subset(df, svid%in% svids)$ac))
  expect_true(all(rowSums(ac)==ac.tot))
  expect_true(!is.null(colnames(ac)))
  expect_true(!is.null(rownames(ac)))
  expect_gt(sum(ac), 0)
})

test_that("Count alleles errors is arg problems", {
  expect_error({ac = countAlleles('../truthrefalt.vcf', sv.sites=list(s='s'))}, 'not found')
  expect_error({ac = countAlleles('../truth.refalt.vcf', sv.sites='not a list')}, 'not a list')
})
