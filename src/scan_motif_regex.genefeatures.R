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
--motifRegexFile='<path to motif regex file>' A file contaiing motif regexes. Two columns: 1. motif name 2. regex
--geneIdsDir='<path to file containing gene ids'> Ensembl format gene ids - one id per line
--help                              display this help text and exit

Example:
Rscript scan_RBPmotifs.genefeatures.R \
    --gtfFile=~/Desktop/datasets/Ensembl75.hg19.gtf \
    --genomeVersion='hg19' \
    --motifRegexFile=/Users/buyar/Desktop/work/useful-pipelines/ARE_motifs.txt \
    --geneIdsDir=~/Desktop/home/buyar/collaborations/mikuda/motif_discovery/geneSets/Group1_WT_up_EDC4_IKK_up.txt"

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

if(!("motifRegexFile" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide a file containing motif regular expressions\n")
}

if(!("geneIdsDir" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide a directory containing a lists of ENSEMBL gene ids\n")
}

gtfFile = argsL$gtfFile
genomeVersion = argsL$genomeVersion
motifRegexFile = argsL$motifRegexFile
geneIdsDir = argsL$geneIdsDir

if(!genomeVersion %in% c('hg19', 'ce6', 'mm9', 'dm3')){
  cat(helpCommand,"\n")
  stop("Error: The supported genome versions are hg19, ce6, mm9 and dm3\n")
}

# load libraries
library('RCAS')
library('stringr')
#data(gff)
#gtf <- gff

# features: GRanges object 
# seqs: DNAStringSet object with names
printHitsTable <- function (features, seqs, motifRegex, outFileBase) {
  features.df <- data.frame('names' = names(features), 'chr' = seqnames(features), 'st' = start(features), 'end' = end(features), 
                            'strand' = strand(features), stringsAsFactors = FALSE)
  hits <- lapply(motifRegex[,c('regex')], FUN = function(x) { str_count(paste(seqs), x)})
  hits.df <- as.data.frame(hits)
  colnames(hits.df) <- motifRegex$regex
  hits.df$totalHits <- rowSums(hits.df)
  hitsTable <- cbind(features.df, hits.df)
  write.table(x = hitsTable, file = paste0(outFileBase, '.motifHitsTable.tsv'), quote = FALSE, sep='\t')
  names(seqs) <- paste0(names(seqs), '.UTR.', c(1:length(seqs)))
  writeXStringSet(x = seqs, filepath = paste0(outFileBase, '.UTRseqs.fasta'))
}

calculateEnrichedMotifs <- function (motifRegex, features, geneIdsFile, genomeVersion) {

  outFileBase <- gsub(pattern = '.txt', replacement = '', x = basename(geneIdsFile))
  geneIds <- scan(geneIdsFile, what = character())
  ##Get foreground and background features
  reducedFeatures <- GenomicRanges::reduce(GenomicRanges::split(x=features, f = features$gene_id))
  foregroundFeatures <- unlist(reducedFeatures[names(reducedFeatures) %in% geneIds])
  allGenes <- names(reducedFeatures)
  backgroundGenes <- sample(allGenes[!allGenes %in% geneIds], 1000)
  backgroundFeatures <- unlist(reducedFeatures[names(reducedFeatures) %in% backgroundGenes])

  fgSeqs <- RCAS::extractSequences(queryRegions = foregroundFeatures, genomeVersion = genomeVersion)
  names(fgSeqs) <- names(foregroundFeatures)
  
  printHitsTable(features = foregroundFeatures, seqs = fgSeqs, motifRegex = motifRegex, outFileBase = paste0(outFileBase, '.fg'))
  
  fgList <- do.call(c, split(x = fgSeqs, f = names(fgSeqs)))
  fgSeqs <- lapply(fgList, toString)
  
  fgHits <- lapply(motifRegex[,c('regex')], FUN = function(x) { str_count(paste(fgSeqs), x)} )
  
  bgSeqs <-  RCAS::extractSequences(queryRegions = backgroundFeatures, genomeVersion = genomeVersion)
  names(bgSeqs) <- names(backgroundFeatures)
  
  printHitsTable(features = backgroundFeatures, seqs = bgSeqs, motifRegex = motifRegex, outFileBase = paste0(outFileBase, '.bg'))
  
  bgList <- do.call(c, split(x = bgSeqs, f = names(bgSeqs)))
  bgSeqs <- lapply(bgList, toString)
  
  bgHits <- lapply(motifRegex[,c('regex')], FUN = function(x) { str_count(paste(bgSeqs), x)} )
  
  #Normalize number of hits per 1000 bps. 
  fgHitsNormalized <- lapply(fgHits, FUN = function (x) { round(x/as.numeric(paste(lapply(fgSeqs, nchar))) * 1000, 2)})
  bgHitsNormalized <- lapply(bgHits, FUN = function (x) { round(x/as.numeric(paste(lapply(bgSeqs, nchar))) * 1000, 2)})
  
  pvalsGreater <- unlist(lapply(c(1:10), function (x) { 
    resGreater <- wilcox.test(fgHitsNormalized[[x]], bgHitsNormalized[[x]], alternative = 'greater')
    return(resGreater$p.value)
  }))
  
  pvalsLess <- unlist(lapply(c(1:10), function (x) { 
    resLess <- wilcox.test(fgHitsNormalized[[x]], bgHitsNormalized[[x]], alternative = 'less')
    return(resLess$p.value)
  }))
  
  
  fgMeans <- unlist(lapply(fgHits, mean))
  bgMeans <- unlist(lapply(bgHits, mean))
  
  fgMeansNormalized <- unlist(lapply(fgHitsNormalized, mean))
  bgMeansNormalized <- unlist(lapply(bgHitsNormalized, mean))
  
  fgSum <- unlist(lapply(fgHits, sum))
  bgSum <- unlist(lapply(bgHits, sum))
  
  fgGenesWithMotif <- unlist(lapply(fgHits, function(x) {sum(x > 0)}))
  bgGenesWithMotif <- unlist(lapply(bgHits, function(x) {sum(x > 0)}))
  
  results <- data.frame(motif = motifRegex$regex, 
                        #fg.Total.Count = fgSum, 
                        #fg.Average = round(fgMeans,1), 
                        fg.Norm.Average = round(fgMeansNormalized,1),
                        #bg.Total.Count = bgSum, 
                        #bg.Average = round(bgMeans,1),
                        bg.Norm.Average = round(bgMeansNormalized,1),
                        fg.PercentWithMotif = round(fgGenesWithMotif/length(fgSeqs)*100,1),
                        bg.PercentWithMotif = round(bgGenesWithMotif/length(bgSeqs)*100,1),
                        #'wilcox.test.pvalue' = pvals,
                        'wilcox.pvalue.greater' = pvalsGreater,
                        'wilcox.pvalue.less' = pvalsLess)
  
  #test for the difference in the UTR lengths of gene groups
  fg.utrLengths <- as.numeric(paste(lapply(fgSeqs, nchar)))
  bg.utrLengths <- as.numeric(paste(lapply(bgSeqs, nchar)))
  p.utrGreater <- wilcox.test(fg.utrLengths, bg.utrLengths, alternative = 'greater')
  p.utrShorter <- wilcox.test(fg.utrLengths, bg.utrLengths, alternative = 'less')
  
  
  sink(file = paste0(outFileBase, '.motif.stats.txt'))
  cat('#Foreground:\n')
  cat('#Number of genes in foreground:',length(unique(names(foregroundFeatures))),'\n')
  cat('#Average length of UTR region per gene:',round(sum(as.numeric(paste(lapply(fgSeqs, nchar))))/length(fgSeqs)),'\n')
  cat('#\n#Background:\n')
  cat('#Number of genes in background:',length(unique(names(backgroundFeatures))),'\n')
  cat('#Average length of UTR region per gene:',round(sum(as.numeric(paste(lapply(bgSeqs, nchar))))/length(bgSeqs)),'\n')
  cat('#\n#Is the UTR length of foreground genes LONGER than the background? - wilcox pvalue (greater):',p.utrGreater$p.value,'\n')
  cat('#Is the UTR length of foreground genes SHORTER than the background? - wilcox pvalue (less):',p.utrShorter$p.value,'\n')
  print(results)
  sink()
}

cat('Reading genome annotations from gtf file\n')
#get genome annotations
gtf <- RCAS::importGtf(filePath = gtfFile)

cat('Extracting txdb features from the gtf file\n')
#find all gene features' genomic coordinates
txdbFeatures <- RCAS::getTxdbFeaturesFromGff(gff = gtf)
#Assign gene ids to the transcripts in the txdbFeatures
ids <- lapply(X = txdbFeatures, function(x) { m = match(x$tx_name, gtf$transcript_id); gtf[m,]$gene_id})
for(i in 1:length(names(txdbFeatures))) {
  txdbFeatures[[i]]$gene_id <- ids[[i]]
}
##

#Parse motif regex from files 
motifRegex <- read.table(motifRegexFile, header=T, sep='\t', stringsAsFactors = FALSE)

for (geneIdsFile in list.files(path = geneIdsDir, pattern = '.txt$', full.names = TRUE)) {
  cat('Analysing ',geneIdsFile,'\n')
#Find known motifs enriched in 3'UTRs
  #cat('Analysing motifs in 3UTRs\n')
  calculateEnrichedMotifs(features = txdbFeatures$threeUTRs,
                       geneIdsFile = geneIdsFile,
                     genomeVersion = genomeVersion, 
                        motifRegex = motifRegex)
}
