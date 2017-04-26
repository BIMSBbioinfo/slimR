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
#' @param iupredPath Path to iupred executable (default:
#'   '/home/buyar/tools/iupred').
#' @param vepPath Path to variant effect predictor script default:
#'   '/home/buyar/.local/bin/variant_effect_predictor.pl',
#' @param clinvarDataURL URL to download clinvar data default:
#'   'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20170404.vcf.gz',
#' @param workingDirectory Path to folder where the database should be created
#' @param nodeN Number of cores to run the database generation.
#' @param overwrite TRUE/FALSE (default: FALSE) Boolean value to decide if the
#'   previously downloaded or processed files be overwritten? (Applies to all
#'   files except objects saved as RDS)
#' @param updateFasta TRUE/FALSE (default: TRUE) Boolean to decide if fasta
#'   download module should be run to update fasta.RDS
#' @param updateGff TRUE/FALSE (default: TRUE) Boolean to decide if gff download
#'   module should be run to update gff.RDS
#' @param updateIupred TRUE/FALSE (default: TRUE) Boolean to decide if disorder
#'   prediction should be run to update iupred.RDS
#' @param updateSlims TRUE/FALSE (default: TRUE) Boolean to decide if motif
#'   prediction should be run to update slims.RDS
#' @param updateVariants TRUE/FALSE (default: TRUE) Boolean to decide if
#'   variation data from (clinvar and humsavar) should be downloaded and updated
#' @param updateSlimChanges TRUE/FALSE (default: TRUE) Boolean to decide if
#'   motif changes should be calculated to update slimChanges.RDS files for
#'   variants
#' @importFrom parallel makeCluster
#' @importFrom doParallel registerDoParallel
#' @export
createDB <- function(uniprotAccessions,
                     iupredPath = '/home/buyar/tools/iupred',
                     vepPath = '/home/buyar/.local/bin/variant_effect_predictor.pl',
                     clinvarDataURL = 'ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20170404.vcf.gz',
                     motifRegex = slimR::motifRegex,
                     workingDirectory = getwd(),
                     nodeN = 1,
                     overwrite = FALSE,
                     updateFasta = TRUE,
                     updateGff = TRUE,
                     updateIupred = TRUE,
                     updateSlims = TRUE,
                     updateVariants = TRUE,
                     updateSlimChanges = TRUE) {

  if(!dir.exists(workingDirectory)){
    dir.create(workingDirectory)
  }

  setwd(workingDirectory)

  if((updateSlims == TRUE | updateSlimChanges == TRUE)
     & updateFasta == FALSE & !file.exists('./slimDB/fasta.RDS')){
    stop("./slimDB/fasta.RDS file is missing and updateFasta is set to FALSE.
         In this condition, updateSlims or updateSlimChanges cannot work.
         Please set updateFasta to TRUE and re-run.
         Or provide a fasta.RDS file in the folder slimDB")
  }

  if(!dir.exists('slimDB')) {
    dir.create('slimDB')
  }

  if (updateFasta == TRUE) {
    fasta <- downloadUniprotFiles(uniprotAccessions = uniprotAccessions,
                                  outDir = workingDirectory,
                                  format = 'fasta',
                                  overwrite = overwrite,
                                  nodeN = nodeN)
    saveRDS(object = fasta, file = './slimDB/fasta.RDS')
  }

  if(updateGff == TRUE) {
    gff <- downloadUniprotFiles(uniprotAccessions = uniprotAccessions,
                                outDir = workingDirectory,
                                format = 'gff',
                                overwrite = overwrite,
                                nodeN = nodeN)
    saveRDS(object = gff, file = './slimDB/gff.RDS')
    rm(gff)
  }

  if(updateIupred == TRUE){
    fastaFiles <- dir(path = './fasta', pattern = '.fasta$', full.names = TRUE)
    ##TODO only pass fastaFiles which match uniprotAccessions
    iupred <- runIUPred(iupredPath = iupredPath,
                        fastaFiles = fastaFiles,
                        overwrite = overwrite,
                        nodeN = nodeN)
    saveRDS(object = iupred, file = './slimDB/iupred.RDS')
    rm(iupred)
  }

  if(updateSlims == TRUE) {
    cl <- parallel::makeCluster(nodeN)
    doParallel::registerDoParallel(cl)

    fasta <- readRDS('./slimDB/fasta.RDS')

    slims <- foreach (i=1:length(fasta), .inorder = TRUE, .verbose = TRUE) %dopar% {
      sequence <- fasta[[i]][1] #sometimes the fasta item may contain multiple sequences
                                #for the same id. Then use the first sequence.
      if(!is.na(sequence)) {
        result <- slimR::searchSLiMs(sequence = sequence,
                                     motifRegex = slimR::motifRegex)
        if(is.null(result)){
          result <- 'No slims found'
        }
      } else {
        result <- 'Fasta file empty'
      }
      result
    }
    names(slims) <- names(fasta)
    saveRDS(object = slims, file = './slimDB/slims.RDS')
    rm(slims)
    stopCluster(cl)
  }

  if(updateVariants == TRUE) {
    fasta <- readRDS('./slimDB/fasta.RDS')
    getClinVarData(url = clinvarDataURL, overwrite = overwrite, parseDownloadedFile = FALSE)
    vcfFilePath <- file.path(workingDirectory, gsub('.gz$', '', basename(clinvarDataURL)))
    runVEP(vepPATH = vepPath, vcfFilePath = vcfFilePath, overwrite = overwrite)
    vepFilePath <- gsub('.vcf$', '.VEP.tsv', vcfFilePath)
    variants <- combineClinVarWithHumsavar(vcfFilePath = vcfFilePath,
                               vepFilePath = vepFilePath,
                               nodeN = nodeN)
    variants <- validateVariants(df = variants, fasta = fasta, nodeN = nodeN)
    saveRDS(object = variants, file = file.path('./slimDB', 'variants.RDS'))
}

  if(updateSlimChanges == TRUE) {

    fasta <- readRDS('./slimDB/fasta.RDS')
    variants <- readRDS('./slimDB/variants.RDS')
    variants <- variants[uniprotAccession %in% uniprotAccessions & validity == 'valid']
    variants$key <- c(1:nrow(variants))
    ##Find Motif Changes by Mutations start###
    diseaseVars <- unique(variants[humsavarVariant == 'Disease' | clinvarVariant == 'Disease'])
    polymorphisms <- unique(variants[!(variants$key %in% diseaseVars$key) &
                                       (humsavarVariant == 'Polymorphism' |
                                          clinvarVariant == 'Polymorphism')])

    slimChangesDisease <- slimR::findMotifChangesMulti(
      sequences = fasta,
      variants = diseaseVars,
      motifRegex = slimR::motifRegex,
      nodeN = nodeN)

    saveRDS(object = slimChangesDisease,
            file = './slimDB/slimChangesDisease.RDS')
    rm(slimChangesDisease)

    slimChangesPolymorphisms <- slimR::findMotifChangesMulti(
      sequences = fasta,
      variants = polymorphisms,
      motifRegex = slimR::motifRegex,
      nodeN = nodeN)
    saveRDS(object = slimChangesPolymorphisms,
            file = './slimDB/slimChangesPolymorphisms.RDS')
    rm(slimChangesPolymorphisms)
  }
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

    downloadFlag <- 0
    if (!file.exists(fileOut) | overwrite == TRUE | (file.exists(fileOut) & file.size(fileOut) == 0)) {
      downloadFlag <- download.file(url = fileUrl, destfile = fileOut, quiet = TRUE, mode = 'w', method = 'wget')
    }

    result <- ''
    if(downloadFlag == 0 & file.size(fileOut) > 0) {
      if(format == 'fasta') {
          result <- paste(Biostrings::readAAStringSet(filepath = fileOut,
                                                      format = 'fasta'))[1]
        } else if (format == 'gff') {
          result <- rtracklayer::import.gff(con = fileOut)
        } else if (format == 'txt') {
          result <- readLines(con = fileOut)
        }
    }
    result
  }
  names(results) <- uniprotAccessions

  emptyResults <- names(results[is.na(results)])
  if(length(emptyResults) > 0) {
    warning(length(emptyResults), " of the downloaded files are empty.
            \tExcluding the following uniprot accessions from the analysis:\n",
            paste0(emptyResults, collapse = '\t'),"\n")
    results <- results[!is.na(results)]
  }
  stopCluster(cl)
  return(results)
}


#' findModifiedSlims
#'
#'
#' Find out if a given slim pattern is only functional if it is
#' post-translationally modified
#'
#' Modified residues are denoted within round brackets according to the
#' convention employed by the ELM annotations e.g. The regex "S..([ST]).."
#' suggests that [ST] residue needs to be modified for the function of this
#' slim. However, round brackets are also used as part of universal regexes to
#' group items. E.g.: motifRegex$TRG_ER_diArg_1 =>
#' "([LIVMFYWPR]R[^YFWDE]{0,1}R)|(R[^YFWDE]{0,1}R[LIVMFYWPR])' So, we need a
#' function to distinguish the usual usage of round brackets from those that
#' suggest a modified residue.
#'
#' @param regex A regular expression pattern from ELM database
#' @export
findModifiedSlims <- function(regex) {
  #find all sub patterns that have round brackets around them
  matches <- slimR::locateAllRegex(sequence = regex, pattern = '\\(.*?\\)')$match
  #loop over each match to find out if there is anything excep Serine, Threonine, or Tyrosines
  result <- c(regex, FALSE, '')
  for (m in matches) {
    #count number of residue positions within each match
    # if the match corresponds to a single residue, then it is a modification site
    variablePositions <- slimR::locateAllRegex(m, '\\[.*?\\]')$match #positions that can be described by square brackets e.g. [IVL]  I or V or L is accepted
    rm <- gsub(pattern = '\\[.*?\\]', '', m)
    singleAApositions <- slimR::locateAllRegex(rm, '[A-Z]')$match
    rm <- gsub(pattern = '[A-Z]', '', rm)
    rm <- gsub(pattern = '[\\(\\)]', '', rm)
    totalSize <- sum(length(variablePositions) + length(singleAApositions))

    if(totalSize == 1 & nchar(rm) == 0) {
      result = c(regex, TRUE, c(variablePositions, singleAApositions))
      break
    }
  }
  return(result)
}
