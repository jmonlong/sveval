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


test_that("Ranges in a region in ggplot", {
  ranges.l = list(
    calls=readSVvcf('../calls.s0.vcf'),
    truth=readSVvcf('../truth.refalt.vcf')
  )
  reg.gr = ranges.l$calls[10]
  pdf('temp.pdf')
  print(plot_ranges(ranges.l, reg.gr))
  dev.off()
  expect_true(file.remove('temp.pdf'))
})

test_that("Ranges in a region in ggplot with SV ids", {
  ranges.l = list(
    calls=readSVvcf('../calls.s0.vcf', keep.ids=TRUE),
    truth=readSVvcf('../truth.refalt.vcf', keep.ids=TRUE)
  )
  reg.gr = ranges.l$calls[10]
  pdf('temp.pdf')
  print(plot_ranges(ranges.l, reg.gr))
  dev.off()
  expect_true(file.remove('temp.pdf'))
})

test_that("Subsect evaluation results by region", {
  bed = data.frame(chr='x', start=c(1e6), end=c(1.49e6))
  reg = GenomicRanges::makeGRangesFromDataFrame(bed)
  eval = subset_eval(res, regions.gr=reg)
  expect_gt(nrow(eval$eval), 0)
  expect_true(any(eval$eval$TP>0))
  expect_true(sum(eval$curve$TP>0)>0)  
})

test_that("Subsect evaluation results by filter", {
  eval = subset_eval(res, accepted.filters='PASS')
  expect_gt(nrow(eval$eval), 0)
  expect_true(any(eval$eval$TP>0))
  expect_true(sum(eval$curve$TP>0)>0)  
})

