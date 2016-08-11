#find gain or loss of SLiMs via point mutations in trans-membrane and non-transmembrane proteins
library('slimR')
#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 3) {
  args <- c("--help")
}

helpCommand = "
download_uniprot.R: given a list of UniProt Accessions, download fasta and gff files for each protein

Arguments:
--unilist=<path to txt file>            List of uniprot accessions; one entry per line
--uniprotDataDir=<path to folder>           Target directory where downloaded uniprot fasta/gff files should be stored
--polymorphicVariantsFile<path to file>    The path to the file that contains human variation data from Uniprot (see help('parseUniprotHumanVariationData') function)
--help                              display this help text and exit

Example:
Rscript mutantSLiMAnalysis \
--unilist=./uniprotList.txt \
--uniprotDataDir=/data/akalin/buyar/collaborations/selbach/data/uniprot/fasta \
--polymorphicVariantsFile=/data/akalin/buyar/collaborations/selbach/data/uniprot/homo_sapiens_variation.txt.gz "


## Help section
if("--help" %in% args) {
  cat(helpCommand, "\n")
  q(save="no")
}

## Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1

if(!("unilist" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide the path to text file containing uniprot accessions\n")
}

if(!("uniprotDataDir" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide a target directory where downloaded uniprot fasta/gff/txt files should be stored\n")
}

if(!("polymorphicVariantsFile" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("Missing polymorphism data\n")
}

unilist <- argsL$unilist
uniprotDataDir <- argsL$uniprotDataDir
polymorphicVariantsFile <- argsL$polymorphicVariantsFile


if(!dir.exists(uniprotDataDir)) {
  stop(uniprotDataDir, 'Folder does not exist')
}

if (!file.exists(polymorphicVariantsFile)) {
  stop(polymorphicVariantsFile, "File does not exist. See help('parseUniprotHumanVariationData')")
}

# read the list of uniprot accessions
uniList <- scan(file = unilist, what = character())

################# Find cytosolic transmembrane regions of proteins
ptm <- proc.time()
cat('Finding cytosolic transmembrane regions of proteins')
pb <- txtProgressBar(min = 0, max = length(uniList), style = 3)
ctmRegions <- list()
if(length(uniList) > 0) {
  for (i in 1:length(uniList)) {
    setTxtProgressBar(pb, i)
    uniAcc <- uniList[i]
    downloadUniprotFiles(uniAcc, outDir = uniprotDataDir, format = 'gff')
    gffFile <- file.path(uniprotDataDir, 'gff', paste0(uniAcc, '.gff'))

    if (file.exists(gffFile)) {
      features <- rtracklayer::import.gff(gffFile)
      ctm <- features[features$type == 'Topological domain',]
      ctm <- ctm[grep(pattern = 'Cytoplasmic', x = ctm$Note),]
      if(length(ctm) > 0) {
        ctmRegions[[length(ctmRegions)+1]] <- ctm
        names(ctmRegions)[length(ctmRegions)] <- uniAcc
      }
    }
  }
}
close(pb)
proc.time() - ptm


################# findMotifChanges due to mutations for a given list of uniprot accessions

variants <- getHumSavar()

ptm <- proc.time()

downloadUniprotFiles(uniprotAccessions = names(ctmRegions), outDir = uniprotDataDir, format = 'fasta')
fastaFiles <- file.path(uniprotDataDir, 'fasta', paste0(names(ctmRegions), '.fasta'))

motifChangesDisease <- list()
motifChangesPolymorphism <- list()

pb <- txtProgressBar(min = 0, max = length(fastaFiles), style = 3)

for (i in 1:length(fastaFiles)) {
  setTxtProgressBar(pb, i)

  fastaFile <- fastaFiles[i]

  #read each fasta file, get the uniprotAcc from file name
  #use slimR::downloadUniprotFiles to get fasta and gff files into a directory

  uniAcc <- gsub(pattern = '.fasta$', replacement = '', x = basename(fastaFile))

  #cat('analysing',uniAcc,'\n')

  sequence <- paste(Biostrings::readAAStringSet(filepath = fastaFile, format = 'fasta'))

  diseaseVariants <- variants[variants$uniprotAcc == uniAcc & variants$variant == 'Disease',]
  if(nrow(diseaseVariants) > 0) {
    changes <- findMotifChanges(sequence = sequence, variants = diseaseVariants, motifRegex = motifRegex)
    if (!is.null(changes)) {
      motifChangesDisease[[length(motifChangesDisease)+1]] <- changes
      names(motifChangesDisease)[length(motifChangesDisease)] <- uniAcc
    }
  }
  polymorphicVariants <- variants[variants$uniprotAcc == uniAcc & variants$variant == 'Polymorphism',]
  if(nrow(polymorphicVariants) > 0) {
    changes <- findMotifChanges(sequence = sequence, variants = polymorphicVariants, motifRegex = motifRegex)
    if (!is.null(changes)) {
      motifChangesPolymorphism[[length(motifChangesPolymorphism)+1]] <- changes
      names(motifChangesPolymorphism)[length(motifChangesPolymorphism)] <- uniAcc
    }
  }
}

close(pb)
proc.time() - ptm

#################TODO: find motif changes in cytosolic regions of transmembrane proteins

#################TODO: compare frequency of loss/gain of motifs in disease vs polymorphism datasets


