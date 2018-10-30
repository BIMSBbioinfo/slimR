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
#' @return A Granges object containing the coordinates of mutated sites in proteins
#' @export
getHumSavar <- function () {
  variantFile <- file.path(getwd(), 'humsavar.txt')
  if (!file.exists(variantFile)) {
    download.file(url = 'www.uniprot.org/docs/humsavar.txt',
                  destfile = variantFile, method = "curl")
  } else {
    warning("humsavar.txt exists at current folder",getwd(),
            ", a new one won't be downloaded. Remove the existing
            file and re-run the function to update the file")
  }

  if(file.exists(paste0(variantFile, '.RDS'))){
    return(readRDS(paste0(variantFile, '.RDS')))
  } else {
    #skip first 30 lines which don't contain mutation data
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
    saveRDS(object = mut, file = paste0(variantFile, ".RDS"))
    return(mut)
  }
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
getClinVarData <- function(url, overwrite = FALSE) {
  destFile <- basename(url)
  if(overwrite == TRUE) {
    download.file(url = url, destfile = destFile)
  } else if (overwrite == FALSE) {
    if(!file.exists(destFile)) {
      download.file(url = url, destfile = destFile)
    }
  }
  if(file.exists(destFile)) {
    if(overwrite == TRUE) {
      gunzipCommand <- paste('gunzip -f',destFile)
    } else {
      gunzipCommand <- paste('gunzip',destFile)
    }
    system(gunzipCommand)
    destFile <- file.path(getwd(), gsub('.gz$', '', destFile))
    return(destFile)
  } else {
    stop("Couldn't find",destFile,"to parse the results.
         Probably the download didn't work\n")
  }
  }

#' getUniprotData
#'
#' Batch download different types of data from Uniprot using
#' REST services
#'
#' @param format File format. Accepted formats are gff, fasta, txt
#' @param organism Organism code e.g. 9606 for human
#' @param reviewed Only download reviewed uniprot entries (default: TRUE)
#' @return
getUniprotData <- function(format, reviewed = TRUE, organism = 9606, overwrite = FALSE) {
  if(!format %in% c('gff', 'fasta', 'txt')) {
    stop("Error: can only download gff, fasta, or txt format files")
  }

  url <- 'https://www.uniprot.org/uniprot/?query='

  if(reviewed == TRUE) {
    url <- paste0(url, "reviewed:yes+AND+")
  }

  url <- paste0(url, "organism:",organism, "&format=",format)

  downloadDate <- paste(unlist(strsplit(date(), ' '))[c(2,3,5)], collapse = "_")
  outFile <- paste(c("uniprot", organism, downloadDate, format), collapse = '.')

  if(file.exists(outFile) & overwrite == FALSE) {
    stop("Destination file already exists: ",outFile,
         "\n Set 'overwrite' to TRUE to overwrite the existing file")
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
#' @param fastaFiles A vector of file paths each pointing to a fasta file
#'   containing the amino acid sequence of a single protein
#' @param outDir The location where the IUPred result files should be stored.
#' @param overwrite TRUE or FALSE (default). Whether the IUPred results should
#'   be overwritten
#' @param nodeN default: 1. Positive integer for the number of parallel cores to use
#' for downloading and processing the files
#' @export
runIUPred <- function (iupredPath,
                       fastaFiles,
                       outDir = getwd(),
                       overwrite = FALSE,
                       nodeN = 1,
                       returnResultsAsList = TRUE) {

  iupredOutDir <- file.path(outDir, 'iupredResults')

  if(!dir.exists(outDir)) {
    dir.create(outDir)
  }
  if(!dir.exists(iupredOutDir)) {
    dir.create(iupredOutDir)
  }

  cl <- parallel::makeCluster(nodeN)
  doParallel::registerDoParallel(cl)

  results <- foreach (i=1:length(fastaFiles), .inorder = TRUE) %dopar% {
    fastaFile <- fastaFiles[i]
    id <- gsub('.fasta', '', basename(fastaFile))
    iupredOut <- file.path(iupredOutDir, paste0(id, '.iupred'))
    if (!file.exists(iupredOut) | overwrite == TRUE) {
      myCommand <- paste(paste0("export IUPred_PATH=",iupredPath, ";"),
                         file.path(iupredPath, 'iupred'),
                         fastaFile, "short >", iupredOut)
      system(command = myCommand)
    }

    if(file.exists(iupredOut)) {
      df <- read.table(file = iupredOut)
      result <- as.numeric(df$V3)
    } else {
      result <- NA
    }
  }
  names(results) <- gsub(pattern = '.fasta$', replacement = '', x = basename(fastaFiles))
  stopCluster(cl)
  return(results)
}


