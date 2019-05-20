context('Graphs')

res = svevalOl('../calls.s0.vcf', '../truth.symb.vcf')

test_that("PR graphs", {
  pdf('temp.pdf')
  ggp.l = plot_prcurve(res$curve)
  expect_warning({tmp = lapply(ggp.l, print)}, 'Removed')
  dev.off()
  expect_true(file.remove('temp.pdf'))
  pdf('temp.pdf')
  ggp.l = plot_prcurve(list(a=res$curve, b=res$curve))
  expect_warning({tmp = lapply(ggp.l, print)}, 'Removed')
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
