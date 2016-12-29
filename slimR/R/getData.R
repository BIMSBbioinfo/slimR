#' createDB
#'
#' Wrapper function to run extractProteinData in multiple nodes for for multiple
#' inputs and create a single database collecting all available information into
#' a single R object
#'
#' @param uniprotAccessions Vector of uniprot accessions (e.g. c('P04637',
#'   'P11166'))
#' @param workingDirectory Path to folder where the database should be created
#' @param nodeN Number of cores to run the database generation.
#' @param overwrite TRUE/FALSE (default: FALSE) if the database is populated in
#'   an existing folder, should the previously generated data be overwritten? Set
#'   to TRUE to update/overwrite existing files/folders
#' @param saveIndividualRDS TRUE/FALSE (default: FALSE) Boolean value to decide if
#' the extracted data for each protein should be saved individually to separate
#' RDS files.
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @export
createDB <- function(iupredPath = '/Users/buyar/Desktop/work/tools/iupred',
                     uniprotAccessions,
                     workingDirectory = getwd(),
                     nodeN = 1,
                     overwrite = FALSE) {
  cl <- parallel::makeCluster(nodeN)
  doParallel::registerDoParallel(cl)

  database <- foreach (i=1:length(uniprotAccessions),
                       .combine = c,
                       .packages = c('GenomicRanges')) %dopar% {
    slimR::extractProteinData(iupredPath = iupredPath,
                                    uni = uniprotAccessions[i],
                                    workingDirectory = workingDirectory,
                                    overwrite = overwrite,
                                    saveAsRDS = TRUE)
    uniprotAccessions[i]
  }
  stopCluster(cl)
}

#' extractProteinData
#'
#' Given a uniprot accession, extract and process all available sequence
#' features including: 1. Fasta sequences 2. Disorder predictions using IUPred
#' 3. Uniprot feature files downloaded in gff format 4. Uniprot variants
#' including disease-causing and polymorphic substitutions 5. Short linear
#' motifs - all substrings matching all available patterns 6. SLiM Changes - The
#' collection of changes in the proteins' SLiM content when variants at 4) are
#' applied to the sequence
#'
#' @param iupredPath The path to the IUPred Disorder score prediction tool
#'   executable
#' @param uni The uniprot accession of the protein for which to extract data
#'   (e.g. P04637)
#' @param workingDirectory The name of the folder to write downloaded/processed
#'   files and RDS objects
#' @param overwrite TRUE/FALSE (default: FALSE). Boolean to decide if downloaded
#'   files should overwrite existing files
#' @param saveAsRDS TRUE/FALSE (default: FALSE). Boolean to decide if all
#'   extracted data should be saved as an RDS file in the working directory
#'   under folder named 'RDS'
#' @return A list object with variable data types directory
#' @export
extractProteinData <- function (iupredPath,
                      uni,
                      workingDirectory = getwd(),
                      overwrite = FALSE,
                      saveAsRDS = FALSE) {

  setwd(dir = workingDirectory)
  variants <- getHumSavar()

  variants <- variants[GenomicRanges::seqnames(variants) == uni,]

  downloadUniprotFiles(uniprotAccessions = uni,
                       format = 'fasta',
                       overwrite = overwrite)

  downloadUniprotFiles(uniprotAccessions = uni,
                       format = 'gff',
                       overwrite = overwrite)


  db <- list('sequence' = NA, 'iupredScore' = NA, 'diseaseVars' = NA,
             'polymorphicVars' = NA, 'uniprotFeatures' = NA, 'slims' = NA,
             'slimChangesDisease' = NA, 'slimChangesPolymorphism' = NA)

  fastaFile <- file.path(getwd(), 'fasta', paste0(uni, '.fasta'))

  if (file.exists(fastaFile)) {
    sequence <- paste(Biostrings::readAAStringSet(filepath = fastaFile,
                                                  format = 'fasta'))

    if(length(sequence) == 0){
      warning("Fasta file for ",uni,"is empty.
      Stopped collecting other features for this id.
              Check file:",fastaFile)
      return(0)
    } else {
      db$sequence <- sequence
    }

    iupredResults <- runIUPred(iupredPath = iupredPath,
                               fastaFiles = c(fastaFile),
                               overwrite = overwrite)

    if(uni %in% names(iupredResults)) {
      db$iupredScore <- iupredResults[[uni]]
      if (nrow(db$iupredScore) == 0) {
        warning("IUPred file for ",uni,"is empty.")
        db$iupredScore = NA
      }
    }

    slims <- slimR::searchSLiMs(sequence = sequence,
                                motifRegex = slimR::motifRegex)
    if(length(slims) > 0){
      db$slims <- slims
    }

    db$diseaseVars <- variants[variants$variant == 'Disease',]
    db$polymorphicVars <- variants[variants$variant == 'Polymorphism',]

    if(length(db$diseaseVars) > 0){
      db$slimChangesDisease <- slimR::findMotifChanges(
        sequence = sequence,
        variants = unique(
          data.frame('wtAA' = db$diseaseVars$wtAA,
                     'mutAA' = db$diseaseVars$mutAA,
                     'pos' = GenomicRanges::start(db$diseaseVars))),
        motifRegex = slimR::motifRegex)
    }

    if (length(db$polymorphicVars) > 0) {
      db$slimChangesPolymorphism <- slimR::findMotifChanges(
        sequence = sequence,
        variants = unique(
          data.frame('wtAA' = db$polymorphicVars$wtAA,
                     'mutAA' = db$polymorphicVars$mutAA,
                     'pos' = GenomicRanges::start(db$polymorphicVars))),
        motifRegex = slimR::motifRegex)
    }
  }

  gffFile <- file.path(getwd(), 'gff', paste0(uni, '.gff'))
  if (file.exists(gffFile)) {
    db$uniprotFeatures <- rtracklayer::import.gff(con = gffFile)
    if(length(db$uniprotFeatures) == 0){
      warning("GFF file for ",uni,"seems to be empty or could not be imported
              Check file:",gffFile)
      db$uniprotFeatures <- NA
    }
  }


  if(saveAsRDS == TRUE) {

    if(!dir.exists('RDS')) {
      dir.create('RDS')
    }

    rdsFile <- file.path(workingDirectory, 'RDS', paste0(uni, '.RDS'))
    if(file.exists(rdsFile)) {
      if(overwrite == TRUE) {
        saveRDS(object = db, file = rdsFile)
      }
    } else {
      saveRDS(object = db, file = rdsFile)
    }
  }

  return(db)
  }

#' parseMutation
#'
#' Given a vector of mutation substitutions (e.g. "p.His160Arg")
#' -> split "p.His160Arg" into "H 160 R"
#' @export
parseMutation <- function (mutations) {
  mutations <- gsub(pattern = '^p.', replacement = '', x = mutations)
  pos <- stringr::str_match(pattern = '\\d+', string = mutations)
  df <- data.frame(do.call(rbind, stringr::str_split(pattern = '\\d+', string = mutations)))
  df$pos <- pos
  colnames(df) <- c('wtAA', 'mutAA', 'pos')
  df$wtAA <- slimR::aaTable$oneLetterCode[match(df$wtAA, slimR::aaTable$threeLetterCode)]
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
                  destfile = variantFile)
  } else {
      warning("humsavar.txt exists at current folder",getwd(),
              ", a new one won't be downloaded. Remove the existing
              file and re-run the function to update the file")
  }

  if(file.exists(paste0(variantFile, '.RDS'))){
    return(readRDS(paste0(variantFile, '.RDS')))
  } else {
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

    parsedMut <- parseMutation(mutations = mut$change)

    mut <- cbind(mut, parsedMut)

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
#' run_IUPred.R: given a directory of fasta sequences of proteins; calculate
#' per-base iupred disorder scores. IUPred source code can be dowloaded from
#' here: http://iupred.enzim.hu/Downloads.php After unpacking the source code,
#' cd to the src directory. Compile the code cc iupred.c -o iupred
#'
#' @param iupredPath The path to the folder containing the source code of the
#'   IUPred tool
#' @param fastaFiles A vector of file paths each pointing to a fasta file
#'   containing the amino acid sequence of a single protein
#' @param outDir The location where the IUPred result files should be stored.
#' @param overwrite TRUE or FALSE (default). Whether the IUPred results should
#'   be overwritten
#' @param returnResultsAsList TRUE (default) or FALSE. Whether the result file should be
#'   parsed and returned as a list of data.frame objects
#' @export
runIUPred <- function (iupredPath, fastaFiles, outDir = getwd(), overwrite = FALSE,
                       returnResultsAsList = TRUE) {

  iupredOutDir <- file.path(outDir, 'iupredResults')

  iupredResultFiles <- c()

  if(!dir.exists(outDir)) {
    dir.create(outDir)
  }
  if(!dir.exists(iupredOutDir)) {
    dir.create(iupredOutDir)
  }

  for (i in 1:length(fastaFiles)) {
    fastaFile <- fastaFiles[i]

    iupredOut <- file.path(iupredOutDir, paste0(gsub(pattern = '.fasta$',
                                                     replacement = '.iupred',
                                                     x = basename(fastaFile))))
    iupredResultFiles <- c(iupredResultFiles, iupredOut)
    if (!file.exists(iupredOut) | overwrite == TRUE) {
      myCommand <- paste(paste0("export IUPred_PATH=",iupredPath, ";"),
                         file.path(iupredPath, 'iupred'),
                         fastaFile, "long >", iupredOut)
      system(command = myCommand)
    }
  }
  if (returnResultsAsList == TRUE) {
    results <- readIUPred(iupredResultFiles = iupredResultFiles)
    names(results) <- gsub(pattern = '.iupred$', replacement = '', x = names(results))
    return (results)
  }
}

#' readIUpred
#'
#' Given a list of files each containing iupred results for a single protein,
#' return a list of data frames where each data frame holds the iupred results.
#'
#' @param A vector of file paths containing the results of IUPred
#' @return A list of data.frame objects
readIUPred <- function(iupredResultFiles) {
  results <- list()
  for(f in iupredResultFiles) {
    df <- read.table(f)
    colnames(df) <- c('pos', 'AA', 'score')
    results[[length(results)+1]] <- df
    names(results)[length(results)] <- basename(f)
  }
  return(results)
}

#' downloadUniprotFiles
#'
#' Given a uniprot accession number and a file format, download the annotation
#' file for a protein in uniprot in the given format. e.g. For human p53
#' protein, the uniprot accession is P04637. This function can be used to
#' download the fasta sequence of this protein from
#' http://www.uniprot.org/uniprot/P04637.fasta. This function is mostly useful
#' for bulk downloads when there are many uniprot accessions to be downloaded.
#'
#' @param uniprotAccessions A vector of uniprot accession numbers
#' @param outDir Path to folder where downloaded files should be stored.
#' @param format The format of files to be downloaded. Options are fasta, gff,
#'   and txt.
#' @param overwrite TRUE or FALSE (default). Whether the existing files should
#'   be overwritten with a new download.
#'
#' @examples
#' ids <- c('P04637', 'P11166', 'P06400')
#' downloadUniprotFiles(uniprotAccessions = ids, outDir = getwd(),
#'                                 format = 'fasta', overwrite = FALSE)
#'
#' @export
downloadUniprotFiles <- function (uniprotAccessions, outDir = getwd(), format, overwrite = FALSE) {
  if (!format %in% c('fasta', 'gff', 'txt')) {
    stop("Uniprot files can only be downloaded in fasta, gff or txt formats.",
         format,"is not a valid option.")
  }
  subOutDir <- file.path(outDir, format)

  if(!dir.exists(outDir)) {
    dir.create(outDir)
  }
  if(!dir.exists(subOutDir)) {
    dir.create(subOutDir)
  }

  for (i in 1:length(uniprotAccessions)) {
    id <- uniprotAccessions[i]
    fileUrl <- paste0("http://www.uniprot.org/uniprot/", id, ".", format)
    fileOut <- file.path(subOutDir, paste0(id, '.', format))

    if (!file.exists(fileOut) | overwrite == TRUE) {
      download.file(url = fileUrl, destfile = fileOut, quiet = TRUE, mode = 'w')
    }
  }
}

