#given a list of Ensembl gene ids and genome version,
#find out if there are any motifs enriched in certain gene features
#(e.g. in UTRs)

#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 2) {
  args <- c("--help")
}

helpCommand = "
RCAS: run RCAS pipeline for a given ENSEMBL annotation file and BED formatted input file

Arguments:
--gtfFile=annotation.gtf     ENSEMBL GTF file
--genome_version=VERSION            genome version; supported values for VERSION are
'hg19'
--help                              display this help text and exit

Example:
RCAS --gtfFile=/data/akalin/Base/Annotation/hg19/EnsemblGenes/142812_EnsemblFTP_hg19_EnsemblGenes.gtf \
--genome_version=hg19"



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

if(!("genome_version" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide the genome version: choose between hg19, mm9, ce6 or dm3\n")
}

gtfFile = argsL$gtfFile
genomeVersion = argsL$genomeVersion

if(!genomeVersion %in% c('hg19', 'ce6', 'mm9', 'dm3')){
  cat(helpCommand,"\n")
  stop("Error: The supported genome versions are hg19, ce6, mm9 and dm3\n")
}

# load libraries
library('RCAS')

gtf <- RCAS::importGtf(filePath = gtfFile)
txdbFeatures <- RCAS::getTxdbFeaturesFromGff(gff = gtf)


# seqs <- lapply(X = txdbFeatures,
#                FUN = function (x) {
#                  RCAS::extractSequences(queryRegions = x[sample(length(x), size = 10000)],
#                                         genomeVersion = genomeVersion)})



