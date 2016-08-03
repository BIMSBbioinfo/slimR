#1. Collect arguments
args <- commandArgs(TRUE)

## Default setting when no arguments passed
if(length(args) < 2) {
  args <- c("--help")
}

helpCommand = "
run_tmhmm.R: given a directory of fasta sequences of proteins; calculate existence probability of trans-membrane helices

tmhmm source code can be downloaded here: http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?tmhmm

Arguments:
--fastaDir=<path to fasta directory>          Directory containing fasta files, one sequence per file
--outDir=<path to output directory>           Target directory where iupred results for each protein should be stored
--overwrite=<yes/no>             Whether the file should be downloaded even if it already exists (default:no)
--help                              display this help text and exit

Example:
Rscript run_IUPred.R \
--fastaDir=./uniprot/fasta \
--outDir=./tmhmmResults \
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

if(!("fastaDir" %in% argsDF$V1)) {
  cat(helpCommand, "\n")
  stop("provide the path to directory with fasta files\n")
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

fastaDir <- argsL$fastaDir
outDir <- argsL$outDir

iupredOutDir <- file.path(outDir, 'tmhmmResults')

if(!dir.exists(outDir)) {
  dir.create(outDir)
}
if(!dir.exists(iupredOutDir)) {
  dir.create(iupredOutDir)
}

# read the list of uniprot accessions
files <- dir(path = fastaDir, pattern = '.fasta$', full.names = TRUE)
pb <- txtProgressBar(min = 0, max = length(files), style = 3)

for (i in 1:length(files)) {
  setTxtProgressBar(pb, i)
  fastaFile <- files[i]

  iupredOut <- file.path(iupredOutDir, paste0(gsub(pattern = '.fasta$', replacement = '.tmhmm', x = basename(fastaFile))))

  if (!file.exists(iupredOut) | overwrite == 'yes') {
    system(command = paste("tmhmm", fastaFile, ">", iupredOut))
  }
}
close(pb)

