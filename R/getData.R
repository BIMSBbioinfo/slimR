#' parseMutation
#'
#' Given a vector of mutation substitutions (e.g. "p.His160Arg")
#' -> split "p.His160Arg" into "H 160 R"
#' @export
parseMutation <- function (mutations) {
  mutations <- gsub(pattern = '^p.', replacement = '', x = mutations)
  pos <- stringr::str_match(pattern = '\\d+', string = mutations)
  df <- data.frame(do.call(rbind, stringr::str_split(pattern = '\\d+', string = mutations)))
  df$pos <- as.numeric(pos)
  colnames(df) <- c('wtAA', 'mutAA', 'pos')
  df$wtAA <- slimR::aaTable$oneLetterCode[match(df$wtAA, slimR::aaTable$threeLetterCode, nomatch = )]
  df$mutAA <- slimR::aaTable$oneLetterCode[match(df$mutAA, slimR::aaTable$threeLetterCode)]
  return(df)
}


#' getHumSavar
#'
#' Download and parse protein mutation data from UniProt
#' @param outdir Path to the folder, where the downloaded file should be saved.
#' @return A Granges object containing the coordinates of mutated sites in
#'   proteins
#' @export
getHumSavar <- function (outdir) {
  variantFile <- file.path(outdir, 'humsavar.txt')
  if (!file.exists(variantFile)) {
    download.file(url = 'www.uniprot.org/docs/humsavar.txt',
                  destfile = variantFile)
  } else {
    warning("humsavar.txt exists at the target location",variantFile,
            ", a new one won't be downloaded. Remove the existing
            file and re-run the function to update the file")
  }

  #skip first 50 lines which don't contain mutation data
  dat <- readLines(con = variantFile)[-(1:50)]

  #grep the lines with relevant variant data
  mut <- grep(pattern = "Polymorphism|Disease|Unclassified",
              x = dat,
              perl = TRUE,
              value = TRUE)
  mut <- data.frame(
    do.call(rbind,
            stringr::str_split(string = mut,
                               n = 7,
                               pattern = '\\s+')))

  colnames(mut) <- c('geneName', 'uniprotAccession',
                     'FTId', 'change',
                     'variant', 'dbSNP', 'diseaseName')

  parsedMut <- parseMutation(mutations = mut$change)

  mut <- cbind(mut, parsedMut)

  #some of the mutations may not have been parsed correctly (due to errors in
  #humsavar file). So, some amino acids may have NA values after conversion.
  #such entries are althogether deleted in the merged data frmae)
  mut <- mut[!(is.na(mut$wtAA) | is.na(mut$mutAA) | is.na(mut$pos)),]

  mut <- GenomicRanges::makeGRangesFromDataFrame(df = mut,
                                                 keep.extra.columns = TRUE,
                                                 ignore.strand = TRUE,
                                                 seqnames.field = 'uniprotAccession',
                                                 start.field = 'pos', end.field = 'pos')
  return(mut)
}


### Functions to download and map ClinVar data to Uniprot sequences
#' getClinVarData
#'
#' This function will fetch the clinvar data from the given url and parse the
#' contents of the downloaded file.
#'
#' @param url The url to the ftp location of the clinvar dataset (tab delimited
#'   variant summary file) (e.g.
#'   ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)
#'
#' @param overwrite default: FALSE, boolean value to decide if a fresh download
#'   should overwrite the existing file
#' @return Path to the downloaded and unzipped file
#' @export
getClinVarData <- function(url = 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz',
                           overwrite = FALSE) {
  destFile <- basename(url)
  if(overwrite == TRUE) {
    download.file(url = url, destfile = destFile)
  } else if (overwrite == FALSE) {
    if(!file.exists(destFile)) {
      download.file(url = url, destfile = destFile)
    }
  }
  return(destFile)
}

#' getUniprotData
#'
#' Batch download different types of data from Uniprot using
#' REST services
#'
#' @param format File format. Accepted formats are gff, fasta, txt
#' @param organism Organism code e.g. 9606 for human
#' @param reviewed Only download reviewed uniprot entries (default: TRUE)
#' @export
getUniprotData <- function(format, reviewed = TRUE, organism = 9606, overwrite = FALSE) {
  if(!format %in% c('gff', 'fasta', 'txt')) {
    stop("Error: can only download gff, fasta, or txt format files")
  }

  url <- 'https://www.uniprot.org/uniprot/?query='

  if(reviewed == TRUE) {
    url <- paste0(url, "reviewed:yes+AND+")
  }

  url <- paste0(url, "organism:",organism, "&format=",format)

  downloadDate <- paste(unlist(strsplit(date(), ' '))[c(2,4,6)], collapse = '_')
  outFile <- paste(c("uniprot", organism, downloadDate, format), collapse = '.')

  if(file.exists(outFile) & overwrite == FALSE) {
    warning("Destination file already exists: ",outFile,
         "\n Set 'overwrite' to TRUE to overwrite the existing file")
    return(outFile)
  }

  download.file(url = url, destfile = outFile)
  return(outFile)
}

#' getElmClasses
#'
#' Scrape ELM classes table from elm.eu.org
#' @export
getElmClasses <- function () {
  urlData <- RCurl::getURL(url = "http://elm.eu.org/elms")
  htmlData <- XML::htmlParse(urlData)
  tables <- XML::readHTMLTable(doc = htmlData)
  elmClasses <- tables[[1]]
  colnames(elmClasses) <- gsub(pattern = ' ',
                               replacement = '_',
                               x = colnames(elmClasses))

  motifRegex <- as.list(as.vector(elmClasses$RegEx))
  names(motifRegex) <- as.vector(elmClasses$ELM_Identifier)

  return(motifRegex)
}


#' getElmInstances
#'
#' Using REST services, download ELM instances table from elm.eu.org
#' @param query (default: * -all instances-) any query term to search the
#'   database for instances
#' @param instanceLogic (default: true+positive) other options: true+negative,
#'   false+positive, false+negative
#' @param taxon (default: homo+sapiens) other option examples: mus+musculus,
#'   drosophila+melanogaster etc.
#'
#' @export
getElmInstances <- function (query = '*',
                             instanceLogic = 'true+positive',
                             taxon = 'homo+sapiens') {
  outFile <- file.path(getwd(), 'elmInstances.tsv')
  if(!file.exists(outFile)) {
    myUrl <- paste0('http://elm.eu.org/instances.tsv?q=',query,
                    '&instance_logic=',instanceLogic,
                    '&taxon=',taxon
                    )
    download.file(url = myUrl, destfile = outFile)
  } else {
    warning("elmInstances.tsv exists at current folder",getwd(),
            ", a new one won't be downloaded. Remove the existing
            file and re-run the function to update the file")
  }
  df <- read.table(outFile, sep = '\t', skip = 5, header = TRUE)
  elmInstances <- GenomicRanges::makeGRangesFromDataFrame(
    df = df,
    keep.extra.columns = TRUE,
    ignore.strand = TRUE,
    seqnames.field = 'Primary_Acc',
    start.field = 'Start',
    end.field = 'End')
  return(elmInstances)
}

#' runIUPred
#'
#' run IUPred tool to get disorder propensity scores of a given protein sequence
#' in fasta format.
#'
#' Given a fasta file containing one or more protein sequences; calculate
#' per-base iupred disorder scores. IUPred source code can be dowloaded from
#' here: http://iupred.enzim.hu/Downloads.php After unpacking the source code,
#' cd to the src directory. Compile the code with "cc iupred.c -o iupred"
#'
#' @param iupredPath The path to the folder containing the source code of the
#'   IUPred tool
#' @param fastaFile A fasta file
#'   containing the amino acid sequences
#' @param nodeN default: 1. Positive integer for the number of parallel cores to use
#' for downloading and processing the files
#' @export
runIUPred <- function (iupredPath,
                       fastaFile,
                       nodeN = 1) {

  fasta <- Biostrings::readAAStringSet(filepath = fastaFile, format = 'fasta')
  cl <- parallel::makeCluster(nodeN)
  parallel::clusterExport(cl = cl, varlist = c('fasta', 'iupredPath'), envir = environment())
  results <- pbapply::pblapply(cl = cl, X = 1:length(fasta), FUN = function(i) {
    require(Biostrings)
    sequence <- fasta[i]
    #write sequence to a temp file
    id <- paste0(c('tmp_', sample(1:9, 10, replace = T)),  collapse = '')
    tmpFasta <- file.path(getwd(), paste0(id, '.fasta'))
    Biostrings::writeXStringSet(x = sequence, filepath = tmpFasta, format = 'fasta')

    #run iupred
    tmpOut <- file.path(getwd(), paste0(id, '.iupred'))
    myCommand <- paste(paste0("export IUPred_PATH=",iupredPath, ";"),
                       file.path(iupredPath, 'iupred'),
                       tmpFasta, "short >", tmpOut)

    system(command = myCommand)

    if(file.exists(tmpOut)) {
      df <- read.table(file = tmpOut)
      result <- as.numeric(df$V3)
    } else {
      result <- NA
    }
    #clean up tmp files
    file.remove(tmpFasta)
    file.remove(tmpOut)
    return(result)
  })
  stopCluster(cl)
  names(results) <- names(fasta)
  return(results)
}


#' getPFAM
#'
#' Download pfam domain annotations of complete proteomes
#'
#' @param organism Organism taxonomy id (e.g. human: 9606)
#' @return A Granges object containing the coordinates PFAM domains in protein
#'   sequences
#' @export
getPFAM <- function(organism = 9606, pfam_version = "Pfam30.0") {

  url <- paste0('ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/',
                pfam_version,
                '/proteomes/',
                organism,
                ".tsv.gz")

  outFile <- file.path(getwd(), paste0(pfam_version, '.', organism, '.tsv.gz'))

  if(file.exists(outFile)) {
    stop("PFAM annotation data already exists at:",outFile,"\n")
  }

  download.file(url = url, destfile = outFile)

  pfam <- data.table::fread(outFile, sep = '\t', header = F)
  colnames(pfam) <- c('seqname', 'start', 'end', 'env_start', 'env_end',
                      'pfam_acc', 'pfam_name', 'pfam_type',
                      'hmm_start', 'hmm_end', 'hmm_length', 'bit_score',
                      'e_value', 'clan'
                      )
  pfam <- GenomicRanges::makeGRangesFromDataFrame(pfam,
                                                  keep.extra.columns = TRUE)

  return(pfam)
}


