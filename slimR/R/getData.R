#' createDB
#'
#' Given a vectore of uniProt accessions, extract and process all available
#' sequence features and save the resulting objects in folders and RDS objects
#' for: 1. Fasta sequences 2. Disorder predictions using IUPred 3. Uniprot
#' feature files downloaded in gff format 4. Uniprot variants including
#' disease-causing and polymorphic substitutions 5. Short linear motifs - all
#' substrings matching all available patterns 6. SLiM Changes due to disease
#' causing mutations - The collection of changes in the proteins' SLiM content
#' when diesase-causing variants at 4) are applied to the sequence 7. SLiM
#' Changes due to polymorphisms - The collection of changes in the proteins'
#' SLiM content when polymorphisms at 4) are applied to the sequence
#'
#' @param uniprotAccessions Vector of uniprot accessions (e.g. c('P04637',
#'   'P11166'))
#' @param workingDirectory Path to folder where the database should be created
#' @param nodeN Number of cores to run the database generation.
#' @param overwrite TRUE/FALSE (default: FALSE) Boolean value to decide if the
#'   previously downloaded or processed files be overwritten? (Applies to all
#'   files except objects saved as RDS)
#' @param updateDB TRUE/FALSE (default: TRUE) Boolean to decide if pre-existing
#'   RDS files generated with this function should be updated/overwritten at
#'   each run.
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @export
createDB <- function(uniprotAccessions,
                     iupredPath = '/home/buyar/tools/iupred',
                     motifRegex = slimR::motifRegex,
                     workingDirectory = getwd(),
                     nodeN = 1,
                     overwrite = FALSE,
                     updateDB = TRUE) {

  if(!dir.exists(workingDirectory)){
    dir.create(workingDirectory)
  }

  setwd(workingDirectory)

  if(!dir.exists('slimDB')) {
    dir.create('slimDB')
  }

  fasta <- downloadUniprotFiles(uniprotAccessions = uniprotAccessions,
                                outDir = workingDirectory,
                                format = 'fasta',
                                overwrite = overwrite,
                                nodeN = nodeN)
  if(updateDB == TRUE) {
    saveRDS(object = fasta, file = './slimDB/fasta.RDS')
  }

  gff <- downloadUniprotFiles(uniprotAccessions = uniprotAccessions,
                              outDir = workingDirectory,
                              format = 'gff',
                              overwrite = overwrite,
                              nodeN = nodeN)
  #rm(gff)
  if(updateDB == TRUE) {
    saveRDS(object = gff, file = './slimDB/gff.RDS')
  }

  ##Disorder prediction start##
  fastaFiles <- dir(path = './fasta', pattern = '.fasta$', full.names = TRUE)
  iupred <- runIUPred(iupredPath = '~/tools/iupred',
                      fastaFiles = fastaFiles,
                      overwrite = overwrite,
                      nodeN = nodeN)
  #  rm(iupred)
  if(updateDB == TRUE) {
    saveRDS(object = iupred, file = './slimDB/iupred.RDS')
  }
  ##Disorder prediction end##

  # The rest of the datasets are not necessary to run
  # if updateDB is not set to TRUE
  if(updateDB == TRUE) {
    ##Motif scanning start#####
    cl <- parallel::makeCluster(nodeN)
    doParallel::registerDoParallel(cl)

    slims <- foreach (i=1:length(fasta), .inorder = TRUE) %dopar% {
      result <- slimR::searchSLiMs(sequence = fasta[[i]],
                                   motifRegex = slimR::motifRegex)
    }
    names(slims) <- names(fasta)
    saveRDS(object = slims, file = './slimDB/slims.RDS')
    rm(slims)
    stopCluster(cl)
    ##Motif scanning start#####


    ##Find Motif Changes by Mutations start###
    variants <- subset(x = as.data.frame(getHumSavar()),
                       select = c('seqnames', 'start',
                                  'wtAA', 'mutAA', 'variant'))
    colnames(variants) <- c('uniprotAccession', 'pos',
                            'wtAA', 'mutAA', 'variant')

    diseaseVars <- unique(variants[variants$variant == 'Disease',])
    polymorphisms <- unique(variants[variants$variant == 'Polymorphism',])

    ##disease mutations###

    slimChangesDisease <- slimR::findMotifChangesMulti(
      sequences = fasta,
      variants = diseaseVars,
      motifRegex = slimR::motifRegex,
      nodeN = nodeN)

    saveRDS(object = slimChangesDisease, file = './slimDB/slimChangesDisease.RDS')
    rm(slimChangesDisease)

    ## polymorphisms##
    slimChangesPolymorphisms <- slimR::findMotifChangesMulti(sequences = fasta,
                                                             variants = polymorphisms,
                                                             motifRegex = slimR::motifRegex,
                                                             nodeN = nodeN)
    saveRDS(object = slimChangesPolymorphisms, file = './slimDB/slimChangesPolymorphisms.RDS')
    #rm(slimChangesPolymorphisms)
    ##Find Motif Changes by Mutations end###
  }
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
  df$pos <- as.numeric(pos)
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
#' Given a directory of fasta sequences of proteins; calculate
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
#' @param nodeN default: 1. Positive integer for the number of parallel cores to use
#' for downloading and processing the files
#' @examples
#' ids <- c('P04637', 'P11166', 'P06400')
#' downloadUniprotFiles(uniprotAccessions = ids, outDir = getwd(),
#'                                 format = 'fasta', overwrite = FALSE)
#'
#' @export
downloadUniprotFiles <- function (uniprotAccessions,
                                  outDir = getwd(),
                                  format,
                                  overwrite = FALSE,
                                  nodeN = 1) {
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

  cl <- parallel::makeCluster(nodeN)
  doParallel::registerDoParallel(cl)

  results <- foreach (i=1:length(uniprotAccessions), .inorder = TRUE) %dopar% {
    id <- uniprotAccessions[i]
    fileUrl <- paste0("http://www.uniprot.org/uniprot/", id, ".", format)
    fileOut <- file.path(subOutDir, paste0(id, '.', format))

    if (!file.exists(fileOut) | overwrite == TRUE) {
      download.file(url = fileUrl, destfile = fileOut, quiet = TRUE, mode = 'w')
    }

    if(file.exists(fileOut)) {
      if(format == 'fasta') {
        result <- paste(Biostrings::readAAStringSet(filepath = fileOut,
                                                      format = 'fasta'))
      } else if (format == 'gff') {
        result <- rtracklayer::import.gff(con = fileOut)
      } else if (format == 'txt') {
        result <- readLines(con = fileOut)
      } else {
        result <- NA
      }
    } else {
      result <- NA
    }
    result
  }
  names(results) <- uniprotAccessions
  stopCluster(cl)
  return(results)
}

