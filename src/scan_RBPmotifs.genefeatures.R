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
--genome_version=VERSION            genome version; supported values for VERSION are
'hg19'
--geneIdsFile='gene_ids.txt'        A txt file containing one ENSEMBL gene id per line
--knownMotifsDir='<path to RBPmap dir>' Directory containing pssm matrices for known RBP motifs (/home/buyar/datasets/known-motifs/RBP_PSSMs)
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

if(!("geneIdsFile" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide a file containing a list of ENSEMBL gene ids\n")
}

gtfFile = argsL$gtfFile
genomeVersion = argsL$genomeVersion
geneIdsFile = argsL$geneIdsFile

if(!genomeVersion %in% c('hg19', 'ce6', 'mm9', 'dm3')){
  cat(helpCommand,"\n")
  stop("Error: The supported genome versions are hg19, ce6, mm9 and dm3\n")
}

# load libraries
library('RCAS')
data(gff)
gtf <- gff
#gtf <- RCAS::importGtf(filePath = gtfFile)

parseMotifPWMs <- function (motifsPath) {
  files <- list.files(path = motifsPath, pattern = '.txt$', full.names = T)
  pwms <- list()
  for (f in files) {
    motif <- gsub(pattern = '.txt$', replacement = '', x = basename(f))
    mm <- read.table(f, row.names = 1, header=T)
    tmm <- t(mm)
    #replace uracil with T
    rownames(tmm) <- gsub(pattern = 'U', replacement = 'T', x = rownames(tmm))
    pwms[[motif]] <- tmm
  }
  return(pwms)
}


geneIds <- scan(file = geneIdsFile, what = character())

#txdbFeatures <- lapply(X = gtfSubsetList, FUN = function (x) { RCAS::getTxdbFeaturesFromGff(gff = x) })
txdbFeatures <- RCAS::getTxdbFeaturesFromGff(gff = gtf)

#Assign gene ids to the transcripts in the txdbFeatures
ids <- lapply(X = txdbFeatures, function(x) { m = match(x$tx_name, gtf$transcript_id); gtf[m,]$gene_id})
for(i in range(length(names(txdbFeatures)))) {
  txdbFeatures[[i]]$gene_id <- ids[[i]]
}
##

features <- txdbFeatures$threeUTRs
reducedFeatures <- GenomicRanges::reduce(GenomicRanges::split(x=features, f = features$gene_id))

#features for the input gene list
foregroundFeatures <- reducedFeatures[names(reducedFeatures) %in% geneIds]

#features for the randomly selected list of genes (equal to the number of genes in the foreground list)
allGenes <- names(reducedFeatures)
backgroundGenes <- sample(x = allGenes[!allGenes %in% geneIds], size = 1000)

backgroundFeatures <- reducedFeatures[backgroundGenes]
bgSeqs <- RCAS::extractSequences(queryRegions = unlist(backgroundFeatures), genomeVersion = genomeVersion)
fgSeqs <- RCAS::extractSequences(queryRegions = unlist(foregroundFeatures), genomeVersion = genomeVersion)

fgHits <- sum(lengths(lapply(X=fgSeqs, FUN = function(x) {Biostrings::matchPWM(pwm = tmm, subject = x, min.score = '90%')})) > 0)
bgHits <- sum(lengths(lapply(X=bgSeqs, FUN = function(x) {Biostrings::matchPWM(pwm = tmm, subject = x, min.score = '90%')})) > 0)


