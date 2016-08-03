#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 2) {
  args <- c("--help")
}

helpCommand = "
download_uniprot.R: given a list of UniProt Accessions, download fasta and gff files for each protein

Arguments:
--unilist=unilist.txt            List of uniprot accessions; one entry per line
--outDir=./uniprotData           Target directory where downloaded files should be stored
--overwrite=<yes/no>             Whether the file should be downloaded even if it already exists (default:no)
--help                              display this help text and exit

Example:
Rscript download_uniprot.R \
--unilist=./uniprotList.txt \
--outDir=./uniprotData \
--overwrite=no"

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

if(!("outDir" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide a target output directory\n")
}

if(!("overwrite" %in% argsDF$V1)) {
  overwrite <- 'no'
} else {
  overwrite <- argsL$overwrite
}

unilist <- argsL$unilist
outDir <- argsL$outDir

fastaOutDir <- file.path(outDir, 'fasta')
gffOutDir <- file.path(outDir, 'gff')

if(!dir.exists(outDir)) {
  dir.create(outDir)
}
if(!dir.exists(fastaOutDir)) {
  dir.create(fastaOutDir)
}
if (!dir.exists(gffOutDir)) {
  dir.create(gffOutDir)
}

# read the list of uniprot accessions
uni <- scan(file = unilist, what = character())
pb <- txtProgressBar(min = 0, max = length(uni), style = 3)

for (i in 1:length(uni)) {
  setTxtProgressBar(pb, i)
  id <- uni[i]
  gffUrl <- paste0("http://www.uniprot.org/uniprot/", id, ".gff")
  gffOut <- file.path(gffOutDir, paste0(id, '.gff'))

  fastaUrl <- paste0("http://www.uniprot.org/uniprot/", id, ".fasta")
  fastaOut <- file.path(fastaOutDir, paste0(id, '.fasta'))

  if (!file.exists(gffOut) | overwrite == 'yes') {
    download.file(url = gffUrl, destfile = gffOut, quiet = TRUE, mode = 'w')
  }

  if (!file.exists(fastaOut) | overwrite == 'yes') {
    download.file(url = fastaUrl, destfile = fastaOut, quiet = TRUE, mode = 'w')
  }
}
close(pb)





