context("Testing data retrieval functions\n")
library(slimR)


hs <- getHumSavar(outdir = getwd())
test_that("humsavar variants", {
          expect_equal(length(hs) > 78000, TRUE)
          expect_equal(colnames(S4Vectors::mcols(hs)), c('geneName', 'FTId', 'change', 'variant',
                                              'dbSNP', 'diseaseName', 'wtAA', 'mutAA'))
          expect_is(hs, 'GRanges')
  })


cv <- getClinVarData()
#read first 100 lines
cv_head <- data.table::fread(text = readLines(cv, n = 100))
test_that("clinvar_variants", {
          expect_equal(cv, "variant_summary.txt.gz")
          expect_equal(file.exists(cv), TRUE)
          expect_equal(ncol(cv_head), 31)
          expect_equal(colnames(cv_head)[1], "#AlleleID")
  })

#test validateVariants functions
glutFastaFile <- system.file("extdata", "glut.fasta", package = 'slimR')
glutFasta <- Biostrings::readAAStringSet(glutFastaFile)
data("glutMutations")
N <- 10
val <- validateVariants(df = glutMutations[1:N,], fasta = glutFasta, nodeN = 1)
test_that("validate_variants", {
  expect_equal(nrow(val), N)
  expect_equal(ncol(val), ncol(glutMutations) + 2)
  expect_equal(val$dbSNP, glutMutations[1:N,]$dbSNP)
})



#download Tobacco Rattle Virus proteome (it has only 6 genes)
#https://www.uniprot.org/proteomes/UP000001669
up <- getUniprotData(outDir = getwd(), format = 'fasta',
                     reviewed = TRUE, organism = 652939)
up_fasta <- Biostrings::readAAStringSet(up)

test_that("get uniprot data - fasta", {
          expect_equal(length(up_fasta), 6)
          expect_is(up_fasta, "AAStringSet")
          expect_equal(sum(grepl('MVP_', names(up_fasta))), 1)
          })

up_gff_file <- getUniprotData(outDir = getwd(), format = 'gff',
               reviewed = TRUE, organism = 652939)
test_that("get uniprot data - gff", {
          expect_equal(file.exists(up_gff_file), TRUE)
  })

test_that("get uniprot data format", {
          expect_error(getUniprotData(format = 'foo'))
  })

elm_classes <- getElmClasses()
elm2pfam <- getELMdomainInteractions()

test_that("get ELM data", {
          expect_is(elm_classes, 'list')
          expect_equal(elm_classes$LIG_FHA_1, "..(T)..[ILV].")
          expect_is(elm2pfam, 'data.frame')
          expect_equal(colnames(elm2pfam), c('ELM_identifier', 'Interaction_Domain_Id',
                                             'Interaction_Domain_Name'))
  })

pfam <- getPFAM(organism = 9606, pfam_version = 'Pfam30.0')
pfamClans <- droplevels(getPFAMClans())
test_that("get PFAM data", {
          expect_is(pfam, 'GRanges')
          expect_equal(colnames(S4Vectors::mcols(pfam)[1:2]), c('env_start', 'env_end'))
          expect_equal(as.vector(pfamClans[pfamClans$ID == 'TRAF',]$Accession), 'CL0389')
  })



