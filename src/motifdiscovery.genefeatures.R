#given a list of Ensembl gene ids and genome version,
#find out if there are any motifs enriched in certain gene features
#(e.g. in UTRs)

#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 3) {
  args <- c("--help")
}

helpCommand = "
RCAS: run RCAS pipeline for a given ENSEMBL annotation file and BED formatted input file

Arguments:
--gtfFile=annotation.gtf            ENSEMBL GTF file
--genomeVersion=VERSION            genome version; supported values for VERSION are
'hg19'
--geneIdsDir='<path to directory containing gene sets'> One or more files with Ensembl format gene ids - one id per line
--help                              display this help text and exit

Example:
Rscript scan_RBPmotifs.genefeatures.R \
--gtfFile=~/Desktop/datasets/Ensembl75.hg19.gtf \
--genomeVersion='hg19' \
--geneIdsDir=~/Desktop/home/buyar/collaborations/mikuda/motif_discovery/geneSets/"

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

if(!("gtfFile" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide the path to ENSEMBL gtf file\n")
}

if(!("genomeVersion" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide the genome version: choose between hg19, mm9, ce6 or dm3\n")
}

if(!("geneIdsDir" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide a directory containing one ore more files each containing a list of ENSEMBL gene ids\n")
}

gtfFile = argsL$gtfFile
genomeVersion = argsL$genomeVersion
geneIdsDir = argsL$geneIdsDir

if(!genomeVersion %in% c('hg19', 'ce6', 'mm9', 'dm3')){
  cat(helpCommand,"\n")
  stop("Error: The supported genome versions are hg19, ce6, mm9 and dm3\n")
}

# load libraries
library('RCAS')
library('GenomicRanges')
#data(gff)
#gtf <- gff

getReducedFeatureSequences <- function (features, geneIds, genomeVersion) {
  #split features by gene id and reduce to get a non-overlapping list of genomic regions
  #for the corresponding feature type.
  reducedFeatures <- GenomicRanges::reduce(GenomicRanges::split(x=features, f = features$gene_id))
  #get a subset of reduced features
  reducedFeaturesSubset <- unlist(reducedFeatures[names(reducedFeatures) %in% geneIds])
  #remove short sequences
  reducedFeaturesSubset <- reducedFeaturesSubset[GenomicRanges::width(reducedFeaturesSubset) > 8]

  #extract the genomic sequences of the extracted features
  sequences <- RCAS::extractSequences(queryRegions = reducedFeaturesSubset,
                                      genomeVersion = genomeVersion)
  #rename the sequences
  names(sequences) <- paste(names(reducedFeaturesSubset), names(sequences), sep='_')
  return(sequences)
}

cat('Reading genome annotations from gtf file\n')
#get genome annotations
gtf <- RCAS::importGtf(filePath = gtfFile)

myFiles <- list.files(path = geneIdsDir,
                      pattern = '.txt$',
                      full.names = TRUE)

for (geneIdsFile in myFiles) {
  cat('Reading gene ids from',geneIdsFile,'\n')
  ##subset txdbFeature by gene ids from input file
  geneIds <- scan(geneIdsFile, what = character())

  cat('Extracting txdb features from the gtf file\n')
  #find all gene features' genomic coordinates
  txdbFeatures <- RCAS::getTxdbFeaturesFromGff(gff = gtf[gtf$gene_id %in% geneIds,])
  #Assign gene ids to the transcripts in the txdbFeatures
  ids <- lapply(X = txdbFeatures, function(x) { m = match(x$tx_name, gtf$transcript_id); gtf[m,]$gene_id})
  for(i in 1:length(names(txdbFeatures))) {
    txdbFeatures[[i]]$gene_id <- ids[[i]]
  }
  ##
  cat('Extracting 3UTR sequences\n')
  threeUTRsequences <- getReducedFeatureSequences( features = txdbFeatures$threeUTRs,
                                                   geneIds = geneIds,
                                                   genomeVersion = genomeVersion
  )
  outFile <- paste0(gsub(pattern = '.txt', replacement = '.threeUTRs.fasta', basename(geneIdsFile)))
  cat('Writing 3UTR sequences into',outFile,'\n\n')
  Biostrings::writeXStringSet(x = threeUTRsequences, filepath = outFile)

  cat('Extracting 5UTR sequences\n')
  fiveUTRsequences <- getReducedFeatureSequences( features = txdbFeatures$fiveUTRs,
                                                  geneIds = geneIds,
                                                  genomeVersion = genomeVersion
  )
  outFile <- paste0(gsub(pattern = '.txt', replacement = '.fiveUTRs.fasta', basename(geneIdsFile)))
  cat('Writing 5UTR sequences',outFile,'\n\n')
  Biostrings::writeXStringSet(x = fiveUTRsequences, filepath = outFile)
}









