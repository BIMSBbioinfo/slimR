context("Testing data retrieval functions")
library(slimR)
hs <- getHumSavar(outdir = getwd())

test_that("humsavar variants",
          expect_equal(length(hs) > 78000, TRUE),
          expect_equal(colnames(mcols(hs)), c('geneName', 'FTId', 'change', 'variant',
                                              'dbSNP', 'diseaseName', 'wtAA', 'mutAA')),
          expect_is(hs, 'GRanges'))


