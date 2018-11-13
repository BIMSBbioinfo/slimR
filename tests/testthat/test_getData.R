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

#download Tobacco Rattle Virus proteome (it has only 6 genes)
#https://www.uniprot.org/proteomes/UP000001669
up <- getUniprotData(outDir = getwd(), format = 'fasta',
                     reviewed = TRUE, organism = 652939)
up_fasta <- Biostrings::readAAStringSet(up)

test_that("get uniprot data - fasta",
          expect_equal(length(up_fasta), 6),
          expect_is(up_fasta, "AAStringSet"),
          expect_equal(grep('MVP_', names(up_fasta)), 1)
          )

up_gff_file <- getUniprotData(outDir = getwd(), format = 'gff',
               reviewed = TRUE, organism = 652939)
test_that("get uniprot data - gff",
          expect_equal(file.exists(up_gff_file), TRUE))

test_that("get uniprot data format",
          expect_error(getUniprotData(format = 'foo')))


