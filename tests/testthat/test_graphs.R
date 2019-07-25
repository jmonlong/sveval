context('Graphs')

res = svevalOl('../calls.s0.vcf', '../truth.refalt.vcf')

test_that("PR graphs", {
  pdf('temp.pdf')
  ggp.l = plot_prcurve(res$curve)
  dev.off()
  expect_true(file.remove('temp.pdf'))
  pdf('temp.pdf')
  ggp.l = plot_prcurve(list(a=res$curve, b=res$curve))
  dev.off()
  expect_true(file.remove('temp.pdf'))
})

test_that("Per size", {
  pdf('temp.pdf')
  ggp.l = plot_persize(res)
  expect_warning({tmp = lapply(ggp.l, print)}, 'Removed')
  dev.off()
  expect_true(file.remove('temp.pdf'))
  df = plot_persize(res, plot=FALSE)
  expect_gt(nrow(df), 0)
})

test_that("Per region", {
  bed = data.frame(chr='x', start=c(1e5, 7e5), end=c(5e5, 1e6))
  reg = GenomicRanges::makeGRangesFromDataFrame(bed)
  pdf('temp.pdf')
  ggp.l = plot_perregion(res, reg)
  tmp = lapply(ggp.l, print)
  dev.off()
  expect_true(file.remove('temp.pdf'))
  df = plot_perregion(res, reg, plot=FALSE)
  expect_gt(nrow(df), 0)
})
