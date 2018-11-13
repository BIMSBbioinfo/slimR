context("Testing data retrieval functions")
library(slimR)


hs <- getHumSavar(outdir = getwd())
test_that("humsavar variants",
          expect_equal(length(hs) > 78000, TRUE),
          expect_equal(colnames(mcols(hs)), c('geneName', 'FTId', 'change', 'variant',
                                              'dbSNP', 'diseaseName', 'wtAA', 'mutAA')),
          expect_is(hs, 'GRanges'))


cv <- getClinVarData()
#read first 100 lines
cv_head <- data.table::fread(text = readLines(cv, n = 100))
test_that("clinvar_variants",
          expect_equal(cv, "variant_summary.txt.gz"),
          expect_equal(file.exists(cv), TRUE),
          expect_equal(ncol(cv_head), 31),
          expect_equal(colnames(cv_head)[1], "#AlleleID"))





