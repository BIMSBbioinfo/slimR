context("Motif analysis")
library(slimR)

glutFastaFile <- system.file("extdata", "glut.fasta", package = 'slimR')
glutFasta <- Biostrings::readAAStringSet(glutFastaFile)
data("glutMutations")

#get disorder scores for each variant position
glutMutations$iupredScore <- glutIUPred[['P11166']][glutMutations$pos]
#get motif changes for variants in disordered regions
motifChanges <- findMotifChangesMulti(sequences = glutFasta,
                                 variants = glutMutations[glutMutations$iupredScore > 0.4,],
                                 motifRegex = slimR::motifRegex,
                                 nodeN = 1)
test_that("look for motif changes in glut1", {
          expect_equal(nrow(motifChanges), 5)
          expect_equal(ncol(motifChanges), 10)
          expect_equal(table(motifChanges$change)[['gained']], 3)
  })
