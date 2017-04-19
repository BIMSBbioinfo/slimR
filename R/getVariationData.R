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

    colnames(mut) <- c('geneName', 'uniprotAcc',
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
                                                   seqnames.field = 'uniprotAcc',
                                                   start.field = 'pos', end.field = 'pos')
    saveRDS(object = mut, file = paste0(variantFile, ".RDS"))
    return(mut)
  }
  }

#' parseUniprotHumanVariationData
#'
#' Parse the human variation data (homo_sapiens_variation.txt.gz) from Uniprot.
#' The variation annotation contains polymorphism data from COSMIC, 1000GP etc.
#' and mapped to Uniprot proteins.
#' See ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/variants/homo_sapiens_variation.txt.gz
#' TODO: add a parameter to overwrite the rds file generated
#' @param filePath The path to the human variation data file (e.g homo_sapiens_variation.txt.gz) downloaded from Uniprot
#' @return A data.table format table of mutation data
#' @export
parseUniprotHumanVariationData <- function (filePath, outFile = 'parseUniprotHumanVariationData.tsv',
                                            keepColumns = c(2,3,5,6,7,8,14),
                                            sourceDB = '1000Genomes',
                                            consequenceType = 'missense variant') {
  if (!file.exists(paste0(filePath, '.rds'))) {
    con <- file(filePath, 'r')
    out <- file(outFile, 'w')
    while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
      fields <- unlist(strsplit(oneLine, '\t'))
      if(length(fields) == 14) {
        if (fields[5] == consequenceType & grepl(sourceDB, fields[14]) == TRUE) {
          writeLines(text = paste(fields[keepColumns], collapse = '\t'), con = out)
        }
      }
    }
    close(con = con)
    close(con = out)
    dt <- data.table::fread(outFile, sep = '\t', header = FALSE)
    saveRDS(object = dt, file = paste0(filePath, '.rds'))
  } else {
    dt <- readRDS(paste0(filePath, '.rds'))
  }
  return(dt)
}

### Functions to download and map ClinVar data to Uniprot sequences
#' getClinVarData
#'
#' This function will fetch the clinvar data from the given url and parse the contents
#' of the downloaded file.
#'
#' @param url The url to the ftp location of the clinvar dataset
#'  (e.g. ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20170404.vcf.gz)
#' @return A tibble object extracted from the downloaded vcf file
#' @importFrom vcfR read.vcfR
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
    destFile <- gsub('.gz$', '', destFile)
    clinvarData <- vcfR::extract_info_tidy(vcfR::read.vcfR(destFile))
    return(clinvarData)
  } else {
    stop("Couldn't find",destFile,"to parse the results. Probably the download didn't work\n")
  }
}



