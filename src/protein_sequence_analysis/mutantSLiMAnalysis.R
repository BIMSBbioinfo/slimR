#find gain or loss of SLiMs via point mutations in trans-membrane and non-transmembrane proteins
library('slimR')
library('data.table')
library('GenomicRanges')
#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 2) {
  args <- c("--help")
}

helpCommand = "
download_uniprot.R: given a list of UniProt Accessions, download fasta and gff files for each protein

Arguments:
--unilist=<path to txt file>            List of uniprot accessions; one entry per line
--uniprotDataDir=<path to folder>           Target directory where downloaded uniprot fasta/gff files should be stored
--help                              display this help text and exit

Example:
Rscript mutantSLiMAnalysis \
--unilist=./uniprotList.txt \
--uniprotDataDir=/data/akalin/buyar/collaborations/selbach/data/uniprot/fasta "


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

unilist <- argsL$unilist
uniprotDataDir <- argsL$uniprotDataDir


if(!dir.exists(uniprotDataDir)) {
  stop(uniprotDataDir, 'Folder does not exist')
}


# variants A Granges object of variants parsed from Humsavar using getHumSavar() function or a subset of it with the same structure
findMotifChangesBulk <- function (variants) {

  ptm <- proc.time()

  motifChanges <- list()

  #load regular expressions for motifs
  data("motifRegex")

  uniprotAccessions <- unique(as.vector(seqnames(variants)))
  downloadUniprotFiles(uniprotAccessions = uniprotAccessions, outDir = uniprotDataDir, format = 'fasta')

  pb <- txtProgressBar(min = 0, max = length(uniprotAccessions), style = 3)

  for (i in 1:length(uniprotAccessions)) {
    setTxtProgressBar(pb, i)
    uniAcc <- uniprotAccessions[i]
    fastaFile <- file.path(uniprotDataDir, 'fasta', paste0(uniAcc, '.fasta'))

    if(file.exists(fastaFile)) {
      sequence <- paste(Biostrings::readAAStringSet(filepath = fastaFile, format = 'fasta'))

      myVariants <- variants[seqnames(variants) == uniAcc]

      df <- data.frame('wtAA' = as.vector(myVariants$wtAA),
                       'mutAA' = as.vector(myVariants$mutAA),
                       'pos' = start(myVariants))

      changes <- findMotifChanges(sequence = sequence, variants = df, motifRegex = motifRegex)
      if (!is.null(changes)) {
        changes$uniprotAcc <- uniAcc
        motifChanges[[length(motifChanges)+1]] <- changes
        names(motifChanges)[length(motifChanges)] <- uniAcc
      }
    }
  }

  close(pb)
  proc.time() - ptm
  return(motifChanges)
}

# read the list of uniprot accessions
uniList <- scan(file = unilist, what = character())

cat('Finding cytosolic transmembrane regions of proteins','\n')
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
      if (sum(c('type', 'Note') %in% colnames(mcols(features))) == 2) {
        ctm <- features[features$type == 'Topological domain' & grepl(pattern = 'Cytoplasmic', x = features$Note),]
        if(length(ctm) > 0) {
          mcols(ctm) <- subset(x = mcols(ctm), select = c('type', 'Note'))
          ctmRegions[[length(ctmRegions)+1]] <- ctm
          names(ctmRegions)[length(ctmRegions)] <- uniAcc
        }
      }
    }
  }
}
ctmRegions <- unlist(GenomicRanges::GRangesList(ctmRegions))
close(pb)
proc.time() - ptm

#################find variants in cytosolic regions of transmembrane proteins '
cat('finding variants in cytosolic regions of transmembrane proteins\n')
variants <- getHumSavar()
overlaps <- GenomicRanges::findOverlaps(query = variants, subject = ctmRegions)
variantsFiltered <- variants[queryHits(overlaps),]
################# findMotifChanges due to mutations for a given list of uniprot accessions
#motif changes due to disease mutations (that overlap cytosolic regions of trans-membrane proteins)
cat('Finding motif changes due to mutations\n')

diseaseVariants <- variantsFiltered[variantsFiltered$variant == 'Disease',]
mcD <- data.table(do.call(rbind, findMotifChangesBulk(variants = diseaseVariants)))
#motif changes due to polymorphisms  (that overlap cytosolic regions of trans-membrane proteins)
polymorphismVariants <- variantsFiltered[variantsFiltered$variant == 'Polymorphism',]
mcP <- data.table(do.call(rbind, findMotifChangesBulk(variants = polymorphismVariants)))

#################TODO: compare frequency of loss/gain of motifs in disease vs polymorphism datasets
lostD <- mcD[mcD$change == 'lost',]
gainedD <- mcD[mcD$change == 'gained',]

lostP <- mcP[mcP$change == 'lost',]
gainedP <- mcP[mcP$change == 'gained',]

lost <- merge(lostD[,length(change), by=SLiM], lostP[,length(change), by=SLiM], by = 'SLiM', all = TRUE)
colnames(lost) <- c('SLiM', 'diseaseSuccess', 'polymorphismSuccess')
lost[is.na(lost)] <- 0
lost$diseaseFail <- length(diseaseVariants) - lost$diseaseSuccess
lost$polymorphismFail <- length(polymorphismVariants) - lost$polymorphismSuccess

gained <- merge(gainedD[,length(change), by=SLiM], gainedP[,length(change), by=SLiM], by = 'SLiM', all = TRUE)
colnames(gained) <- c('SLiM', 'diseaseSuccess', 'polymorphismSuccess')
gained[is.na(gained)] <- 0
gained$diseaseFail <- length(diseaseVariants) - gained$diseaseSuccess
gained$polymorphismFail <- length(polymorphismVariants) - gained$polymorphismSuccess

lost$pVal <-  apply(lost, 1, FUN = function(x) {
  contingencyTable <- matrix(as.integer(x[2:5]), nrow = 2,
                             dimnames = list(c("disease", "polymorphism"),
                                             c("success", "fail")))
  pVal <- fisher.test(x = contingencyTable, alternative = 'greater')$p.value
})

gained$pVal <- apply(gained, 1, FUN = function(x) {
  contingencyTable <- matrix(as.integer(x[2:5]), nrow = 2,
                             dimnames = list(c("disease", "polymorphism"),
                                             c("success", "fail")))
  pVal <- fisher.test(x = contingencyTable, alternative = 'greater')$p.value
})

lost <- lost[order(pVal)]
gained <- gained[order(pVal)]

write.table(x = lost, file = 'lost_motifs.tsv', quote = FALSE, sep = '\t', row.names = FALSE)
write.table(x = gained, file = 'gained_motifs.tsv', quote = FALSE, sep = '\t', row.names = FALSE)


