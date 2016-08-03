#' getHumSavar
#'
#' Download and parse protein mutation data from UniProt
#' @export
getHumSavar <- function () {
  variantFile <- file.path(getwd(), 'humsavar.txt')
  if (!file.exists(variantFile)) {
    download.file(url = 'www.uniprot.org/docs/humsavar.txt',
                  destfile = variantFile)
  }
  #skip first 30 lines which don't contain mutation data
  dat <- readLines(con = variantFile)[-(1:30)]

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
  colnames(mut) <- c('geneName', 'uniprotAcc',
                     'FTId', 'change',
                     'variant', 'dbSNP', 'diseaseName')
  return(mut)
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
  return(elmClasses)
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
  if(!exists(outFile)) {
    myUrl <- paste0('http://elm.eu.org/instances.tsv?q=',query,
                    '&instance_logic=',instanceLogic,
                    '&taxon=',taxon
                    )
    download.file(url = myUrl, destfile = outFile)
  }
  df <- read.table(outFile, sep = '\t', skip = 5, header = TRUE)
  elmInstances <- GenomicRanges::makeGRangesFromDataFrame(df = df,
                                                keep.extra.columns = TRUE,
                                                ignore.strand = TRUE,
                                                seqnames.field = 'Primary_Acc',
                                                start.field = 'Start',
                                                end.field = 'End')
  return(elmInstances)
}



