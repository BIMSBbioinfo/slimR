#find gain or loss of SLiMs via point mutations in trans-membrane and non-transmembrane proteins
library('slimR')
#setwd('~/Desktop/data/akalin/buyar/collaborations/selbach/mutantSLiMAnalysis')
uniprotFastaDir <- '../data/uniprot/fasta'
#uniprotFastaDir <- '../data/uniprot/sample_fasta/'
uniprotGffDir <- '../data/uniprot/gff'
iupredResultDir <- '../data/iupredResults/'
iupredPath <- '~/Desktop/work/tools/iupred/'
variants <- getHumSavar()
data("motifRegex")



################# findMotifChanges due to mutations for a given list of fasta files in a given directory
if (file.exists(file.path(getwd(), 'motifChanges.rds'))) {
  motifChanges <- readRDS('motifChanges.rds')
} else {
  ptm <- proc.time()

  fastaFiles <- dir(uniprotFastaDir, full.names = TRUE)

  motifChanges <- list()

  pb <- txtProgressBar(min = 0, max = length(fastaFiles), style = 3)

  for (i in 1:length(fastaFiles)) {
    setTxtProgressBar(pb, i)

    fastaFile <- fastaFiles[i]

    #read each fasta file, get the uniprotAcc from file name
    #use slimR::downloadUniprotFiles to get fasta and gff files into a directory

    uniAcc <- gsub(pattern = '.fasta$', replacement = '', x = basename(fastaFile))

    #cat('analysing',uniAcc,'\n')

    sequence <- paste(Biostrings::readAAStringSet(filepath = fastaFile, format = 'fasta'))
    mutations <- variants[variants$uniprotAcc == uniAcc & variants$variant == 'Disease',]
    if(nrow(mutations) > 0) {
      changes <- findMotifChanges(sequence = sequence, variants = mutations, motifRegex = motifRegex)
      if (!is.null(changes)) {
        motifChanges[[length(motifChanges)+1]] <- changes
        names(motifChanges)[length(motifChanges)] <- uniAcc
      }
    }
  }

  close(pb)
  proc.time() - ptm

  saveRDS(object = motifChanges, file = 'motifChanges.rds')
}

#################

################# Find cytosolic transmembrane regions of proteins
if (file.exists(file.path(getwd(), 'ctmRegions.rds'))) {
  ctmRegions <- readRDS('ctmRegions.rds')
} else {
  ptm <- proc.time()
  gffFiles <- dir(path = uniprotGffDir, full.names = TRUE)
  pb <- txtProgressBar(min = 0, max = length(gffFiles), style = 3)
  ctmRegions <- list()
  if(length(gffFiles) > 0) {
    for (i in 1:length(gffFiles)) {
      setTxtProgressBar(pb, i)
      gffFile <- gffFiles[i]
      uniAcc <- gsub(pattern = '.gff$', replacement = '', x = basename(gffFile))
      features <- rtracklayer::import.gff(gffFile)
      ctm <- features[features$type == 'Topological domain',]
      ctm <- ctm[grep(pattern = 'Cytoplasmic', x = ctm$Note),]
      if(length(ctm) > 0) {
        ctmRegions[[length(ctmRegions)+1]] <- ctm
        names(ctmRegions)[length(ctmRegions)] <- uniAcc
      }
    }
  }
  close(pb)
  proc.time() - ptm
  saveRDS(object = ctmRegions, file = 'ctmRegions.rds')
}


################# Find disorder scores of proteins
if (file.exists(file.path(getwd(), 'iupredScores.rds'))) {
  iupredScores <- readRDS(file = 'iupredScores.rds')
} else {
  iupredScores <- runIUPred(iupredPath = iupredPath, outDir = iupredResultDir, fastaFiles = dir(path = uniprotFastaDir, full.names = TRUE))
  saveRDS(object = iupredScores, file = 'iupredScores.rds')
}

polymorphisms <- parseUniprotHumanVariationData(filePath = '../data/uniprot/homo_sapiens_variation.txt.gz')

