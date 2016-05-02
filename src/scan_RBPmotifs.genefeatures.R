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
--knownMotifsDir='<path to RBPmap dir>' Directory containing pssm matrices for known RBP motifs (/home/buyar/datasets/known-motifs/RBP_PSSMs)
--geneIdsFile='<path to file containing gene ids'> Ensembl format gene ids - one id per line
--help                              display this help text and exit

Example:
Rscript scan_RBPmotifs.genefeatures.R \
    --gtfFile=~/Desktop/datasets/Ensembl75.hg19.gtf \
    --genomeVersion='hg19' \
    --knownMotifsDir=~/Desktop/home/buyar/datasets/known-motifs/RBP_PSSMs \
    --geneIdsFile=~/Desktop/home/buyar/collaborations/mikuda/motif_discovery/geneSets/Group1_WT_up_EDC4_IKK_up.txt"

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

if(!("knownMotifsDir" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide a file containing PSSM values for known RBP motifs\n")
}

if(!("geneIdsFile" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide a file containing a list of ENSEMBL gene ids\n")
}

gtfFile = argsL$gtfFile
genomeVersion = argsL$genomeVersion
knownMotifsDir = argsL$knownMotifsDir
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

#features: GRanges object
#PWMs pwms returned from parseMotifPWMs
countMotifHits <- function (PWMs, features, genomeVersion) {

  cat('extracting sequences\n')
  Seqs <- RCAS::extractSequences(queryRegions = features, genomeVersion = genomeVersion)

  pwmHits <- c()
  for (i in 1:length(PWMs)) {
    tmm <- PWMs[[i]]
    motif <- names(PWMs)[i]
    cat('scanning for motif',motif)
    nhits <- sum(lengths(lapply(X=Seqs, FUN = function(x) {Biostrings::matchPWM(pwm = tmm, subject = x)})) > 0)
    cat('=> found ',nhits,'hits\n')
    pwmHits <- c(pwmHits, nhits)
  }
  return (pwmHits)
}

txdbFeatures <- RCAS::getTxdbFeaturesFromGff(gff = gtf)

#Assign gene ids to the transcripts in the txdbFeatures
ids <- lapply(X = txdbFeatures, function(x) { m = match(x$tx_name, gtf$transcript_id); gtf[m,]$gene_id})
for(i in 1:length(names(txdbFeatures))) {
  txdbFeatures[[i]]$gene_id <- ids[[i]]
}
##

##Get PWMs from files
PWMs <- parseMotifPWMs(motifsPath = knownMotifsDir)

calculateEnrichedMotifs <- function (features, geneIds, genomeVersion, PWMs) {
  ##Get foreground and background features
  reducedFeatures <- GenomicRanges::reduce(GenomicRanges::split(x=features, f = features$gene_id))
  foregroundFeatures <- unlist(reducedFeatures[names(reducedFeatures) %in% geneIds])
  allGenes <- names(reducedFeatures)
  backgroundGenes <- sample(allGenes[!allGenes %in% geneIds], 100)
  backgroundFeatures <- unlist(reducedFeatures[names(reducedFeatures) %in% backgroundGenes])

  fgHits <- countMotifHits(PWMs = PWMs, features = foregroundFeatures, genomeVersion = genomeVersion)
  bgHits <- countMotifHits(PWMs = PWMs, features = backgroundFeatures, genomeVersion = genomeVersion)

  fgSize <- length(foregroundFeatures)
  bgSize <- length(backgroundFeatures)

  pvals <- c()
  for (i in 1:length(PWMs)) {
    f <- fgHits[i]
    b <- bgHits[i]
    contingencyTable <- matrix(c(f, b, fgSize - f, bgSize - b), nrow = 2,
                               dimnames = list(c("foreground", "background"),
                                               c("found", "not_found")))
    p  <- stats::fisher.test(contingencyTable, alternative = "greater")$p.value
    pvals <- c(pvals, p)
  }

  results <- data.frame('motifs' = names(PWMs),
                        'fgHits' = fgHits,
                        'fgFraction' = round(fgHits/fgSize,2),
                        'bgHits' = bgHits,
                        'bgFraction' = round(bgHits/bgSize,2),
                        'pValue' = pvals)
  return(results[order(results$pValue),])
}

results <- calculateEnrichedMotifs(features = txdbFeatures$threeUTRs,
                                   geneIds = scan(file = geneIdsFile,
                                                  what = character()),
                                   genomeVersion = genomeVersion, PWMs = PWMs
                                   )
outFile <- paste0(gsub(pattern = '.txt', replacement = '', basename(geneIdsFile)), '.motifscanresults.threeUTRs.tsv')

write.table(x = results, file = outFile, quote = FALSE, row.names = FALSE, sep = '\t')

results <- calculateEnrichedMotifs(features = txdbFeatures$fiveUTRs,
                                   geneIds = scan(file = geneIdsFile,
                                                  what = character()),
                                   genomeVersion = genomeVersion,
                                   PWMs = PWMs)
outFile <- paste0(gsub(pattern = '.txt', replacement = '', basename(geneIdsFile)), '.motifscanresults.fiveUTRs.tsv')

write.table(x = results, file = outFile, quote = FALSE, row.names = FALSE, sep = '\t')


